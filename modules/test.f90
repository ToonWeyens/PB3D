!------------------------------------------------------------------------------!
!   Test some routines and functions                                           !
!------------------------------------------------------------------------------!
module test
    use num_vars, only: dp, max_str_ln, pi
    use output_ops, only: writo, lvl_ud, print_ar_1, print_ar_2, write_out
    use input_ops, only: yes_no
    use time, only: start_time, stop_time
    use str_ops, only: i2str, r2str, r2strt

    implicit none
    private
    public test_repack, test_write_out, test_mesh_cs, test_metric_transf, &
        &test_ang_B, test_ext_var, test_B, test_VMEC_norm_deriv, &
        &test_arr_mult, test_VMEC_conv_FHM, test_calc_RZL, test_calc_T_VF, &
        &test_calc_inv_met, test_det, test_inv

contains
    subroutine test_repack
        ! VMEC variable has structure (1:mnmax, 1:n_r)
        ! output variable should have (1:n_r, 0:mpol-1, -ntor:ntor)
        use fourier_ops, only: repack
        integer :: n_rB, mpolB, ntorB, mnmaxB
        real(dp), allocatable :: xmB(:), xnB(:)
        real(dp), allocatable :: varinB(:,:)
        real(dp), allocatable :: varoutB(:,:,:)
            
        call writo('test repack?')
        if(yes_no(.false.)) then
            call lvl_ud(1)
            n_rB = 2
            mpolB = 2
            ntorB = 3
            mnmaxB = (ntorB+1)+(2*ntorB+1)*(mpolB-1)
            write(*,*) 'mnmaxB = ', mnmaxB
            allocate(xmB(mnmaxB))
            allocate(xnB(mnmaxB))
            allocate(varinB(mnmaxB,n_rB))
            allocate(varoutB(0:mpolB-1,-ntorB:ntorB,1:n_rB))
            write(*,*) 'size(varinB) = ', size(varinB,1), size(varinB,2)
            write(*,*) 'size(varoutB) = ', size(varoutB,1), size(varoutB,2), &
                &size(varoutB,3)
                
            xmB = [0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
            xnB = [0.0, 1.0, 2.0, 3.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0]
            varinB(:,1) = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
            varinB(:,2) = 2*varinB(:,1)
            varoutB = repack(varinB,mnmaxB,n_rB,mpolB,ntorB,xmB,xnB)
                
            !write(*,*) 'varinB(:,1) = ', varinB(:,1)
            !write(*,*) 'varinB(:,2) = ', varinB(:,2)
            write(*,*) 'varoutB(0,:,1) = ', varoutB(0,:,1)
            write(*,*) 'varoutB(1,:,1) = ', varoutB(1,:,1)
            write(*,*) 'varoutB(0,:,2) = ', varoutB(0,:,2)
            write(*,*) 'varoutB(1,:,2) = ', varoutB(1,:,2)
            call lvl_ud(-1)
            
            call writo('Paused... press enter')
            read(*,*)
        end if
    end subroutine

    subroutine test_write_out
        use output_ops, only: write_out
        use num_vars, only: output_i
        use VMEC_vars, only: &
            &mnmax, n_r, rmnc
            
        call writo('test write_out?')
        if(yes_no(.false.)) then
            call writo('Testing write_out')
            call lvl_ud(1)
            call write_out(mnmax, n_r, rmnc, 'rmnc_name', output_i)       
            open (unit=78,file="results.txt",action="write",status="replace")
            call write_out(mnmax, n_r, rmnc, 'rmnc_name', 78)       
            close (78)
            call lvl_ud(-1)
            
            call writo('Paused... press enter')
            read(*,*)
        end if
    end subroutine

    subroutine test_mesh_cs
        use fourier_ops, only: mesh_cs
        use output_ops, only: print_ar_2
     
        real(dp), allocatable :: output(:,:,:)
        real(dp) :: theta, zeta
        integer :: mpol, ntor
        
        allocate(output(0:mpol-1,-ntor:ntor,2))
        theta = 0
        zeta = 0
        
        do
            call writo('test mesh_cs?')
            if(yes_no(.false.)) then
                call writo('Testing write_out')
                call lvl_ud(1)
                
                write(*,'(A)',advance='no') 'mpol, ntor = ' 
                read(*,*) mpol, ntor
                
                write(*,'(A)',advance='no') 'theta, zeta (*pi) = ' 
                read(*,*) theta, zeta
                theta = pi*theta
                zeta = pi*zeta
                
                output =  mesh_cs(mpol,ntor,theta,zeta)
                write(*,*) 'size mesh = ', size(output,2), size(output,3)
                write(*,*) 'mesh_cs(:,:,1) = '
                call print_ar_2(output(:,:,1))
                write(*,*) 'mesh_cs(:,:,2) = '
                call print_ar_2(output(:,:,2))
                call writo('restarting')
                call lvl_ud(-1)
            else
                exit
            end if
            
            call writo('Paused... press enter')
            read(*,*)
        end do
    end subroutine
    
    subroutine test_VMEC_norm_deriv
        use utilities, only: VMEC_norm_deriv
        
        ! local variables
        integer :: loc_max, id, kd
        integer :: n_steps
        real(dp) :: start_step
        real(dp), allocatable :: varin(:)
        real(dp), allocatable :: var1_an(:), var2_an(:), var3_an(:), &
            &var4_an(:), var5_an(:)
        real(dp), allocatable :: var1_nm(:), var2_nm(:), var3_nm(:), &
            &var4_nm(:), var5_nm(:)
        logical :: ind_plots
        real(dp), allocatable :: step_size(:)
        real(dp), allocatable :: max_error(:,:)
        real(dp), allocatable :: mean_error(:,:)
        real(dp), allocatable :: plot_var(:,:)
        integer :: input_type
        
        call writo('test VMEC_norm_deriv?')
        if(yes_no(.false.)) then
            ! read user input
            do
                call writo('how many steps ?')
                read(*,*) n_steps
                if (n_steps.lt.10) then
                    call writo('Choose a value larger than 10')
                    cycle
                else
                    exit
                end if
            end do
            do
                call  writo('starting (max) step size ?')
                read(*,*) start_step
                if (start_step.gt.1E-1_dp) then
                    call writo('Choose a value smaller than 0.1')
                    cycle
                else
                    exit
                end if
            end do
            
            call writo('Individual plots?')
            if(yes_no(.false.)) then
                ind_plots = .true.
            else
                ind_plots = .false.
            end if
            
            do
                call writo('which function?')
                call lvl_ud(1)
                    call writo('1: sin(pi*x) + 0.25*cos(2*pi*x)')
                    call writo('2: 1/800 * (0.5+0.5*x^2-0.005*x^3+0.0.0000001*x^4)')
                    call writo('3: x^7')
                    call writo('4: x^6')
                    call writo('5: x^5')
                    call writo('6: x^4')
                call lvl_ud(-1)
                read(*,*) input_type
                if (input_type.lt.1 .or. input_type.gt.6) then
                    call writo('Choose a value from the range 1-6')
                    cycle
                else
                    exit
                end if
            end do
            
            ! initialize
            allocate(step_size(n_steps)); step_size = 0.0_dp
            allocate(max_error(n_steps,5)); max_error = 0.0_dp
            allocate(mean_error(n_steps,5)); mean_error = 0.0_dp
            allocate(plot_var(2,n_steps)); plot_var = 0.0_dp
            
            call lvl_ud(1)
                
            ! n_steps calculations up to min_step
            do id = 1,n_steps
                loc_max = nint(id/start_step)
                step_size(id) = 2*pi/loc_max
                
                call writo(trim(i2str(id))//'/'//trim(i2str(n_steps))//&
                    &': calculating using '//trim(i2str(loc_max))//' points')
                
                ! set up input function
                if (allocated(varin)) deallocate(varin)
                if (allocated(var1_an)) deallocate(var1_an)
                if (allocated(var2_an)) deallocate(var2_an)
                if (allocated(var3_an)) deallocate(var3_an)
                if (allocated(var4_an)) deallocate(var4_an)
                if (allocated(var5_an)) deallocate(var5_an)
                if (allocated(var1_nm)) deallocate(var1_nm)
                if (allocated(var2_nm)) deallocate(var2_nm)
                if (allocated(var3_nm)) deallocate(var3_nm)
                if (allocated(var4_nm)) deallocate(var4_nm)
                if (allocated(var5_nm)) deallocate(var5_nm)
                
                allocate(varin(loc_max))
                allocate(var1_an(loc_max))
                allocate(var2_an(loc_max))
                allocate(var3_an(loc_max))
                allocate(var4_an(loc_max))
                allocate(var5_an(loc_max))
                allocate(var1_nm(loc_max))
                allocate(var2_nm(loc_max))
                allocate(var3_nm(loc_max))
                allocate(var4_nm(loc_max))
                allocate(var5_nm(loc_max))
                
                select case (input_type)
                    case (1)
                        varin = [(sin((kd-1)*pi/(loc_max-1))&
                            &+0.25*cos((kd-1)*pi/(loc_max-1)),kd=1,loc_max)]
                        var1_an = (0.5)**1*[(cos(pi*(kd-1)/(loc_max-1))&
                            &-0.25*sin(pi*(kd-1)/(loc_max-1)),kd=1,loc_max)]
                        var2_an = (0.5)**2*[(-sin(pi*(kd-1)/(loc_max-1))&
                            &-0.25*cos(pi*(kd-1)/(loc_max-1)),kd=1,loc_max)]
                        var3_an = (0.5)**3*[(-cos(pi*(kd-1)/(loc_max-1))&
                            &+0.25*sin(pi*(kd-1)/(loc_max-1)),kd=1,loc_max)]
                        var4_an = (0.5)**4*[(sin(pi*(kd-1)/(loc_max-1))&
                            &+0.25*cos(pi*(kd-1)/(loc_max-1)),kd=1,loc_max)]
                        var5_an = (0.5)**5*[(cos(pi*(kd-1)/(loc_max-1))&
                            &-0.25*sin(pi*(kd-1)/(loc_max-1)),kd=1,loc_max)]
                    case (2)
                        varin = 1.0_dp/800*[(0.5+0.5*(kd*step_size(id))**2&
                            &-0.005*(kd*step_size(id))**3&
                            &+0.0000001*(kd*step_size(id))**4,kd=1,loc_max)]
                        var1_an = 1.0_dp/800*[((kd*step_size(id))&
                            &-0.015*(kd*step_size(id))**2&
                            &+0.0000004*(kd*step_size(id))**3,kd=1,loc_max)]
                        var2_an = 1.0_dp/800*[(1-0.03*(kd*step_size(id))&
                            &+0.0000012*(kd*step_size(id))**2,kd=1,loc_max)]
                        var3_an = 1.0_dp/800*[(-0.03+0.0000024*(kd*step_size(id))&
                            &,kd=1,loc_max)]
                        var4_an = 1.0_dp/800*[(0.0000024,kd=1,loc_max)]
                        var5_an = 0.0_dp
                    case(3)
                    write(*,*) 'step_size = ', step_size(id)
                        varin = [((kd*step_size(id))**7,kd=1,loc_max)]
                        var1_an = 7*[((kd*step_size(id))**6,kd=1,loc_max)]
                        var2_an = 7*6*[((kd*step_size(id))**5,kd=1,loc_max)]
                        var3_an = 7*6*5*[((kd*step_size(id))**4,kd=1,loc_max)]
                        var4_an = 7*6*5*4*&
                            &[((kd*step_size(id))**3,kd=1,loc_max)]
                        var5_an = 7*6*5*4*3*&
                            &[((kd*step_size(id))**2,kd=1,loc_max)]
                    case(4)
                        varin = [((kd*step_size(id))**6,kd=1,loc_max)]
                        var1_an = 6*[((kd*step_size(id))**5,kd=1,loc_max)]
                        var2_an = 6*5*[((kd*step_size(id))**4,kd=1,loc_max)]
                        var3_an = 6*5*4*[((kd*step_size(id))**3,kd=1,loc_max)]
                        var4_an = 6*5*4*3*&
                            &[((kd*step_size(id))**2,kd=1,loc_max)]
                        var5_an = 6*5*4*3*2*[((kd*step_size(id)),kd=1,loc_max)]
                    case(5)
                        varin = [((kd*step_size(id))**5,kd=1,loc_max)]
                        var1_an = 5*[((kd*step_size(id))**4,kd=1,loc_max)]
                        var2_an = 5*4*[((kd*step_size(id))**3,kd=1,loc_max)]
                        var3_an = 5*4*3*[((kd*step_size(id))**2,kd=1,loc_max)]
                        var4_an = 5*4*3*2*[((kd*step_size(id)),kd=1,loc_max)]
                        var5_an = 5*4*3*2*1
                    case(6)
                        varin = [((kd*step_size(id))**4,kd=1,loc_max)]
                        var1_an = 4*[((kd*step_size(id))**3,kd=1,loc_max)]
                        var2_an = 4*3*[((kd*step_size(id))**2,kd=1,loc_max)]
                        var3_an = 4*3*2*[((kd*step_size(id))**1,kd=1,loc_max)]
                        var4_an = 4*3*2*1
                        var5_an = 0.0_dp
                    case default
                        call writo('ERROR: in test_VMEC_norm_deriv, how did &
                            &you get here?!??!')
                        stop
                end select
                
                ! numerical derivatives
                call VMEC_norm_deriv(varin,var1_nm,1.0_dp/step_size(id),1,1)
                call VMEC_norm_deriv(varin,var2_nm,1.0_dp/step_size(id),2,1)
                call VMEC_norm_deriv(varin,var3_nm,1.0_dp/step_size(id),3,1)
                call VMEC_norm_deriv(varin,var4_nm,1.0_dp/step_size(id),4,1)
                call VMEC_norm_deriv(varin,var5_nm,1.0_dp/step_size(id),5,1)
                
                ! max and mean errors
                max_error(id,1) = maxval(abs(var1_an-var1_nm))/ &
                    &maxval(abs(var1_an))
                max_error(id,2) = maxval(abs(var2_an-var2_nm))/ &
                    &maxval(abs(var2_an))
                max_error(id,3) = maxval(abs(var3_an-var3_nm))/ &
                    &maxval(abs(var3_an))
                max_error(id,4) = maxval(abs(var4_an-var4_nm))/ &
                    &maxval(abs(var4_an))
                max_error(id,5) = maxval(abs(var5_an-var5_nm))/ &
                    &maxval(abs(var5_an))
                
                mean_error(id,1) = sum(abs(var1_an-var1_nm))/loc_max/ &
                    &maxval(abs(var1_an))
                mean_error(id,2) = sum(abs(var2_an-var2_nm))/loc_max/ &
                    &maxval(abs(var2_an))
                mean_error(id,3) = sum(abs(var3_an-var3_nm))/loc_max/ &
                    &maxval(abs(var3_an))
                mean_error(id,4) = sum(abs(var4_an-var4_nm))/loc_max/ &
                    &maxval(abs(var4_an))
                mean_error(id,5) = sum(abs(var5_an-var5_nm))/loc_max/ &
                    &maxval(abs(var5_an))
                
                ! plots
                if (ind_plots) then
                    call write_out(1,loc_max,varin,'input variable')
                    
                    call write_out(1,loc_max,var1_an,'analytical deriv. ord. 1')
                    call write_out(1,loc_max,var1_nm,'numerical deriv. ord. 1')
                    call write_out(1,loc_max,(var1_an-var1_nm)/&
                        &maxval(abs(var1_an)),'diff, deriv. 1')
                    !call write_out(1,loc_max,log10(abs(var1_an-var1_nm)/&
                        !&maxval(abs(var1_an))),'diff, deriv. 1 [log]')
                    
                    call write_out(1,loc_max,var2_an,'analytical deriv. ord. 2')
                    call write_out(1,loc_max,var2_nm,'numerical deriv. ord. 2')
                    call write_out(1,loc_max,(var2_an-var2_nm)/&
                        &maxval(abs(var2_an)),'diff, deriv. 2')
                    !call write_out(1,loc_max,log10(abs(var2_an-var2_nm)/&
                        !&maxval(abs(var2_an))),'diff, deriv. 2 [log]')
                    
                    call write_out(1,loc_max,var3_an,'analytical deriv. ord. 3')
                    call write_out(1,loc_max,var3_nm,'numerical deriv. ord. 3')
                    call write_out(1,loc_max,(var3_an-var3_nm)/&
                        &maxval(abs(var3_an)),'diff, deriv. 3')
                    !call write_out(1,loc_max,log10(abs(var3_an-var3_nm)/&
                        !&maxval(abs(var3_an))),'diff, deriv. 3 [log]')
                    
                    call write_out(1,loc_max,var4_an,'analytical deriv. ord. 4')
                    call write_out(1,loc_max,var4_nm,'numerical deriv. ord. 4')
                    call write_out(1,loc_max,(var4_an-var4_nm)/&
                        &maxval(abs(var4_an)),'diff, deriv. 4')
                    !call write_out(1,loc_max,log10(abs(var4_an-var4_nm)/&
                        !&maxval(abs(var4_an))),'diff, deriv. 4 [log]')
                    
                    call write_out(1,loc_max,var5_an,'analytical deriv. ord. 5')
                    call write_out(1,loc_max,var5_nm,'numerical deriv. ord. 5')
                    call write_out(1,loc_max,(var5_an-var5_nm)/&
                        &maxval(abs(var5_an)),'diff, deriv. 5')
                    !call write_out(1,loc_max,log10(abs(var5_an-var5_nm)/&
                        !&maxval(abs(var5_an))),'diff, deriv. 5 [log]')
                end if
            end do
            
            ! plot errors
            do id = 1,5
                plot_var(1,:) = step_size
                plot_var(2,:) = max_error(:,id)
                call write_out(2,n_steps,plot_var,'f(Delta) = max. error as &
                    &function of step size Delta',comment='for  &
                    &deriv. of order '//trim(i2str(id)))
                plot_var(1,:) = log10(step_size)
                plot_var(2,:) = log10(max_error(:,id))
                call write_out(2,n_steps,plot_var,'f(Delta) = max. error as &
                    &function of step size Delta [log-log]',comment='for  &
                    &deriv. of order '//trim(i2str(id)))
            end do
            do id = 1,5
                plot_var(1,:) = step_size
                plot_var(2,:) = mean_error(:,id)
                call write_out(2,n_steps,plot_var,'f(Delta) = mean error as &
                    &function of step size Delta',comment='for  &
                    &deriv. of order '//trim(i2str(id)))
                plot_var(1,:) = log10(step_size)
                plot_var(2,:) = log10(mean_error(:,id))
                call write_out(2,n_steps,plot_var,'f(Delta) = mean error as &
                    &function of step size Delta [log-log]',comment='for  &
                    &deriv. of order '//trim(i2str(id)))
            end do
            
            call lvl_ud(-1)
        end if
    end subroutine
    
    subroutine test_ext_var
        use utilities, only : ext_var
        
        ! local variables
        integer :: jd, kd
        integer :: n_r
        real(dp), allocatable :: x(:), varin(:), var_num(:),varder(:), &
            &varder_num(:), x_int(:)
        
        call writo('test ext_var')
        if(yes_no(.false.)) then
            call writo('Testing whether extrapolation is working correctly')
            
            n_r = 100
            allocate(x(n_r))
            allocate(varin(n_r))
            allocate(var_num(n_r))
            allocate(varder(n_r))
            allocate(varder_num(n_r))
            allocate(x_int(n_r))
            do kd = 1,n_r
                ! some continuous curve:
                x(kd) = float(kd)*2*pi/n_r
                x_int(kd) = float(2*kd-1)*pi/n_r
                varin(kd) = sin(float(kd)*2*pi/n_r) + &
                    &0.5*cos(3*float(kd)*2*pi/n_r)
                varder(kd) = cos(float(kd)*2*pi/n_r) - &
                    &1.5*sin(3*float(kd)*2*pi/n_r)
            end do
            do kd = 1,n_r-2
                var_num(kd) = ext_var([(varin(jd),jd=kd,kd+2)],&
                    &[(x(jd),jd=kd,kd+2)],x_int(kd),0)
                varder_num(kd) = ext_var([(varin(jd),jd=kd,kd+2)],&
                    &[(x(jd),jd=kd,kd+2)],x(kd),1)
            end do
            call writo('Don''t worry about the last 3 points!!!')
            call write_out(2,n_r,transpose(reshape([x,varin],[n_r,2])),&
                &'function: sin(x)+0.5*cos(3x)')
            call write_out(2,n_r,transpose(reshape([x_int,var_num],[n_r,2])),&
                &'num interp. of function: sin(x)+0.5*cos(3x)')
            call write_out(2,n_r,transpose(reshape([x,varder],[n_r,2])),&
                &'deriv of function: cos(x)-1.5*sin(3x)')
            call write_out(2,n_r,transpose(reshape([x,varder_num],[n_r,2])),&
                &'num interp. of deriv. of function: cos(x)-1.5*sin(3x)')
            call write_out(2,n_r-2,transpose(reshape([x(1:n_r-2),&
                &abs(varder_num(1:n_r-2)-varder(1:n_r-2))],[n_r-2,2])),&
                &'relative error: cos(x)-1.5*sin(3x)')
            
            !! check whether the interpolating polynomial yields the source points
            !do kd = 1,n_r-3
                !var_num(kd) = ext_var([(varin(jd),jd=kd,kd+2)],&
                    !&[(x(jd),jd=kd,kd+2)],x(kd),0)
                !write(*,*) 'kd = ', kd
                !write(*,*) 'var_num = ', var_num(kd)
                !write(*,*) 'varin   = ', varin(kd)
            !end do
        end if
    end subroutine
    
    subroutine test_VMEC_conv_FHM
        use utilities, only: VMEC_conv_FHM
        
        ! local variables
        integer :: id, kd
        integer :: n_max, step, n_r
        integer, parameter :: length = 50
        real(dp), allocatable :: varin(:), varout(:), varoutout(:), vardiff(:)
        real(dp) :: maxerr(length), averr(length), num_points(length)
        real(dp) :: plotvar(2,length)
        logical :: ind_plot, log_plot
        
        call writo('test VMEC_conv_FHM?')
        if(yes_no(.false.)) then
            call writo('Testing whether h2f*f2h = 1')
            call writo('The relative difference between an original FM variable&
                & var and h2f*f2h*var, should be decreasing with increasing &
                &number of radial points')
            call writo('')
            call writo('Do you want the individual plots?')
            ind_plot = yes_no(.false.)
            call writo('')
            do
                call writo('n_max = ?')
                read(*,*) n_max
                if (n_max.lt.4*length) then
                    call writo('n_max has to be larger than or equal to '&
                        &//trim(i2str(4*length))//'...')
                else 
                    exit
                end if
            end do
            call writo('logarithmic plot?')
            log_plot = yes_no(.false.)
            call writo('')
            
            call lvl_ud(1)
            
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
                
                call VMEC_conv_FHM(varin,varout,.true.)
                call VMEC_conv_FHM(varout,varoutout,.false.)
                do kd = 1,n_r
                    vardiff(kd) = 2*(varin(kd)-varoutout(kd))/&
                        &(varin(kd)+varoutout(kd))
                end do
                maxerr(id) = maxval(abs(vardiff))
                averr(id) = sum(abs(vardiff))/size(vardiff)
                
                if(ind_plot) then
                    call write_out(1,n_r,varin,'input with '//trim(i2str(n_r))&
                        &//' radial points ('//trim(i2str(id))//'/'//&
                        &trim(i2str(length))//')')
                    !call write_out(1,n_r,varout,'first output with '//&
                        !&trim(i2str(n_r))//' radial points ('//&
                        !&trim(i2str(id))//'/'//trim(i2str(length))//')')
                    call write_out(1,n_r,varoutout,'output with '//&
                        &trim(i2str(n_r))//' radial points ('//trim(i2str(id))&
                        &//'/'//trim(i2str(length))//')')
                    call write_out(1,n_r,vardiff,'absolute relative difference&
                        & with '//trim(i2str(n_r))//' radial points ('&
                        &//trim(i2str(id))//'/'//trim(i2str(length))//')', &
                        &comment='max = '//trim(r2strt(maxerr(id)))&
                        &//', average: '//trim(r2strt(averr(id))))
                end if
            end do
            
            plotvar(1,:) = num_points
            if (log_plot) then
                plotvar(2,:) = log10(maxerr)
            else
                plotvar(2,:) = maxerr
            end if
            call write_out(2,length,plotvar,'maximum error as a function of &
                &numer of points')
            if (log_plot) then
                plotvar(2,:) = log10(averr)
            else
                plotvar(2,:) = averr
            end if
            call write_out(2,length,plotvar,'average error as a function of &
                &numer of points')
                
            call writo('Paused... press enter')
            read(*,*)
            call lvl_ud(-1)
        end if
    end subroutine
    
    subroutine test_calc_RZL
        use eq_vars, only: calc_RZL, init_eq, eqd_mesh, &
            &VMEC_R, VMEC_Z, theta, theta_H, zeta, zeta_H, n_par
        use VMEC_vars, only: n_r, rmax_surf, rmin_surf, zmax_surf
        
        ! local variables
        integer :: jd
        
        call writo('test calc_RZL?')
        if(yes_no(.false.)) then
            call writo('initializing')
            call lvl_ud(1)
            call init_eq
            call lvl_ud(-1)
            
            call writo('calculations with constant zeta, varying theta')
            call lvl_ud(1)
            if (allocated(zeta)) deallocate(zeta)
            allocate(zeta(n_par,n_r)); zeta = 0.0_dp
            if (allocated(zeta_H)) deallocate(zeta_H)
            allocate(zeta_H(n_par,n_r)); zeta_H = 0.0_dp
            if (allocated(theta)) deallocate(theta)
            allocate(theta(n_par,n_r)); theta = 0.0_dp
            if (allocated(theta_H)) deallocate(theta_H)
            allocate(theta_H(n_par,n_r)); theta_H = 0.0_dp
            
            zeta = 0.4*pi/2
            zeta_H = 0.4*pi/2
            do jd = 1,n_r
                theta(:,jd) = eqd_mesh(n_par, 0.0_dp*pi, 3.0_dp*pi)
                theta_H(:,jd) = eqd_mesh(n_par, 0.0_dp*pi, 3.0_dp*pi)
            end do
            
            call calc_RZL([0,0,0])
            call calc_RZL([1,0,0])
            call calc_RZL([2,0,0])
            call calc_RZL([3,0,0])
            call calc_RZL([0,0,1])
            call calc_RZL([0,0,2])
            call lvl_ud(-1)
            
            call writo('Plotting R')
            call lvl_ud(1)
            call write_out(2,n_par,transpose(reshape([theta(:,5),&
                &VMEC_R(:,5,0,0,0)],[n_par,2])),'R at r=5')
            call write_out(2,n_par,transpose(reshape([theta(:,5),&
                &VMEC_R(:,5,1,0,0)],[n_par,2])),'R_theta at r=5')
            call write_out(2,n_par,transpose(reshape([theta(:,5),&
                &VMEC_R(:,5,2,0,0)],[n_par,2])),'R_theta^2 at r=5')
            call write_out(2,n_par,transpose(reshape([theta(:,5),&
                &VMEC_R(:,5,3,0,0)],[n_par,2])),'R_theta^3 at r=5')
            call write_out(1,n_r,VMEC_R(5,:,0,0,0),'R at theta=5')
            call write_out(1,n_r,VMEC_R(5,:,0,0,1),'R_r at theta=5')
            call write_out(1,n_r,VMEC_R(5,:,0,0,2),'R_r^2 at theta=5')
            call lvl_ud(-1)
            
            call writo('Plotting Z')
            call lvl_ud(1)
            call write_out(2,n_par,transpose(reshape([theta(:,5),&
                &VMEC_Z(:,5,0,0,0)],[n_par,2])),'Z at r=5')
            call write_out(2,n_par,transpose(reshape([theta(:,5),&
                &VMEC_Z(:,5,1,0,0)],[n_par,2])),'Z_theta at r=5')
            call write_out(2,n_par,transpose(reshape([theta(:,5),&
                &VMEC_Z(:,5,2,0,0)],[n_par,2])),'Z_theta^2 at r=5')
            call write_out(2,n_par,transpose(reshape([theta(:,5),&
                &VMEC_Z(:,5,3,0,0)],[n_par,2])),'Z_theta^3 at r=5')
            call write_out(1,n_r,VMEC_Z(5,:,0,0,0),'Z at theta=5')
            call write_out(1,n_r,VMEC_Z(5,:,0,0,1),'Z_r at theta=5')
            call write_out(1,n_r,VMEC_Z(5,:,0,0,2),'Z_r^2 at theta=5')
            call lvl_ud(-1)
            
            ! test whether VMEC given bounds are respected
            if (rmax_surf.gt.0.0_dp) then
                call writo('Checking the bounds of R')
                call lvl_ud(1)
                call within_bounds(VMEC_R(:,:,0,0,0),rmin_surf,rmax_surf)
                call lvl_ud(-1)
            else
                call writo('Not possible to check the bounds of R')
            end if
            if (zmax_surf.gt.0.0_dp) then
                call writo('Checking the bounds of Z')
                call lvl_ud(1)
                call within_bounds(VMEC_Z(:,:,0,0,0),-zmax_surf,zmax_surf)
                call lvl_ud(-1)
            else
                call writo('Not possible to check the bounds of Z')
            end if
        end if
    contains
        ! display whether the calculated table for R or Z fits within the limits
        ! outputted by VMEC
        ! Since the calculations are done following a magnetic field line, it is
        ! actually quite  probable that the  minimum or maximum is  not reached.
        ! Therefore, only hard  warnings are given, when the  minimum or maximum
        ! calculated value is grossly under or above the VMEC minimum or maximum
        subroutine within_bounds(var,min_VMEC,max_VMEC)
            ! input and output
            real(dp), intent(in) :: var(:,:)
            real(dp), intent(in) :: min_VMEC, max_VMEC
            
            ! local variables
            real(dp) :: margin
            real(dp) :: min_frac , max_frac
            
            margin = 5.0E-2_dp                                                  ! 1% for comparison 
            
            ! minimum and maximum of variable
            min_frac = 2*(minval(var)-min_VMEC)/abs(minval(var)+min_VMEC)       ! positive if within bounds
            max_frac = 2*(maxval(var)-max_VMEC)/abs(maxval(var)+max_VMEC)       ! positive if out of bounds
            
            if (min_frac.lt.-margin) then                                       ! too low minimum
                call writo('WARNING: minimum of variable in real angular space &
                    & is lower than VMEC provided minimum by '//&
                    &trim(r2strt(100*min_frac))//'%...')
                call writo(' -> Maybe use an even number of poloidal points, &
                    &or more angular mesh points in general?')
                call writo(' -> Maybe run VMEC with more accuracy?')
            else if (max_frac.gt.margin) then                                   ! too high maximum
                call writo('WARNING: maximum of variable in real angular space &
                    & is greater than VMEC provided maximum by '//&
                    &trim(r2strt(100*max_frac))//'%...')
                call writo(' -> Maybe use more angular mesh points')
                call writo(' -> Maybe run VMEC with more accuracy?')
            else 
                call writo('within bounds')
                return
            end if
        end subroutine
    end subroutine
    
    subroutine test_calc_T_VF
        use eq_ops, only: calc_eq
        use metric_ops, only: T_VF
        use eq_vars, only: q_saf, VMEC_L, flux_p, theta
        use VMEC_vars, only: n_r
        
        ! local variables
        real(dp) :: T_loc(1:n_r,3,3)
        integer :: id, jd
        integer :: case_nr
        integer :: deriv(3)
        
        call writo('test calc_T_VF?')
        if(yes_no(.false.)) then
            do 
                call writo('Which case do you want?')
                call writo('    1: 0th order')
                call writo('    2: 1st order, r derivative')
                call writo('    3: 1st order, theta derivative')
                call writo('    4: 1st order, zeta derivative')
                call writo('    5: 4th order, rthetathetazeta derivative')
                read(*,*) case_nr
                if(case_nr.lt.1 .or. case_nr.gt.5) then
                    call writo('You have to choose from the available options')
                else
                    exit
                end if
            end do
            
            call writo('calculating equilibrium')
            call lvl_ud(1)
            call calc_eq(0.2*pi)           
            call lvl_ud(-1)
            
            call writo('calculating case')
            call lvl_ud(1)
            select case (case_nr)
                case(1)
                    deriv = [0,0,0]
                    call case_1
                case(2)
                    deriv = [1,0,0]
                    call case_2
                case(3)
                    deriv = [0,1,0]
                    call case_3
                case(4)
                    deriv = [0,0,1]
                    call case_4
                case(5)
                    deriv = [1,2,1]
                    call case_5
                case default
                    call writo('ERROR: In test_T_VF, case '//&
                        &trim(i2str(case_nr))//' is not associated with &
                        &anything. How did you get here???')
                    stop
            end select
            call lvl_ud(-1)
            call writo('case calculated')
            
            call writo('plotting case')
            call lvl_ud(1)
            do id = 1,3
                do jd = 1,3
                    call write_out(1,n_r,T_VF(5,:,id,jd,deriv(1),deriv(2),&
                        &deriv(3)),'T_VF('//trim(i2str(id))//','//&
                        &trim(i2str(jd))//')')
                    call write_out(1,n_r,T_loc(:,id,jd),'alternative T_VF('//&
                        &trim(i2str(id))//','//trim(i2str(jd))//')')
                    call write_out(1,n_r,abs(T_loc(:,id,jd)-T_VF(5,:,id,jd,&
                        &deriv(1),deriv(2),deriv(3))),'abs. diff('//&
                        &trim(i2str(id))//','//trim(i2str(jd))//')')
                end do
            end do
            call lvl_ud(-1)
            call writo('case plotted')
        end if
    contains
        subroutine case_1
            T_loc = 0.0_dp
            T_loc(:,1,1) = -q_saf(:,1)*(theta(5,:)+VMEC_L(5,:,0,0,0)) &
                &-q_saf(:,0)*VMEC_L(5,:,1,0,0)
            T_loc(:,1,2) = flux_p(:,1)/(2*pi)
            T_loc(:,1,3) = VMEC_L(5,:,1,0,0)
            T_loc(:,2,1) = -q_saf(:,0)*(1.0_dp+VMEC_L(5,:,0,1,0))
            T_loc(:,2,2) = 0.0_dp
            T_loc(:,2,3) = 1.0_dp + VMEC_L(5,:,0,1,0)
            T_loc(:,3,1) = 1.0_dp - q_saf(:,0)*VMEC_L(5,:,0,0,1)
            T_loc(:,3,2) = 0.0_dp
            T_loc(:,3,3) = VMEC_L(5,:,0,0,1)
        end subroutine
        
        subroutine case_2
            T_loc = 0.0_dp
            T_loc(:,1,1) = -q_saf(:,2)*(theta(5,:)+VMEC_L(5,:,0,0,0)) &
                &-2*q_saf(:,1)*VMEC_L(5,:,1,0,0)&
                &-q_saf(:,0)*VMEC_L(5,:,2,0,0)
            T_loc(:,1,2) = flux_p(:,2)/(2*pi)
            T_loc(:,1,3) = VMEC_L(5,:,2,0,0)
            T_loc(:,2,1) = -q_saf(:,1)*(1.0_dp+VMEC_L(5,:,0,1,0)) &
                &-q_saf(:,0)*VMEC_L(5,:,1,1,0)
            T_loc(:,2,2) = 0.0_dp
            T_loc(:,2,3) = VMEC_L(5,:,1,1,0)
            T_loc(:,3,1) = -q_saf(:,1)*VMEC_L(5,:,0,0,1) &
                &- q_saf(:,0)*VMEC_L(5,:,1,0,1)
            T_loc(:,3,2) = 0.0_dp
            T_loc(:,3,3) = VMEC_L(5,:,1,0,1)
        end subroutine
        
        subroutine case_3
            T_loc = 0.0_dp
            T_loc(:,1,1) = -q_saf(:,1)*(1+VMEC_L(5,:,0,1,0))&
                &-q_saf(:,0)*VMEC_L(5,:,1,1,0)
            T_loc(:,1,2) = 0.0_dp
            T_loc(:,1,3) = VMEC_L(5,:,1,1,0)
            T_loc(:,2,1) = -q_saf(:,0)*VMEC_L(5,:,0,2,0)
            T_loc(:,2,2) = 0.0_dp
            T_loc(:,2,3) = VMEC_L(5,:,0,2,0)
            T_loc(:,3,1) = -q_saf(:,0)*VMEC_L(5,:,0,1,1)
            T_loc(:,3,2) = 0.0_dp
            T_loc(:,3,3) = VMEC_L(5,:,0,1,1)
        end subroutine
        
        subroutine case_4
            T_loc = 0.0_dp
            T_loc(:,1,1) = -q_saf(:,1)*VMEC_L(5,:,0,0,1)&
                &-q_saf(:,0)*VMEC_L(5,:,1,0,1)
            T_loc(:,1,2) = 0.0_dp
            T_loc(:,1,3) = VMEC_L(5,:,1,0,1)
            T_loc(:,2,1) = -q_saf(:,0)*VMEC_L(5,:,0,1,1)
            T_loc(:,2,2) = 0.0_dp
            T_loc(:,2,3) = VMEC_L(5,:,0,1,1)
            T_loc(:,3,1) = -q_saf(:,0)*VMEC_L(5,:,0,0,2)
            T_loc(:,3,2) = 0.0_dp
            T_loc(:,3,3) = VMEC_L(5,:,0,0,2)
        end subroutine
        
        subroutine case_5
            T_loc = 0.0_dp
            T_loc(:,1,1) = -q_saf(:,2)*VMEC_L(5,:,0,2,1)&
                &-2*q_saf(:,1)*VMEC_L(5,:,1,2,1)-q_saf(:,0)*VMEC_L(5,:,2,2,1)
            T_loc(:,1,2) = 0.0_dp
            T_loc(:,1,3) = VMEC_L(5,:,2,2,1)
            T_loc(:,2,1) = -q_saf(:,1)*VMEC_L(5,:,0,3,1)&
                &-q_saf(:,0)*VMEC_L(5,:,1,3,1)
            T_loc(:,2,2) = 0.0_dp
            T_loc(:,2,3) = VMEC_L(5,:,1,3,1)
            T_loc(:,3,1) = -q_saf(:,1)*VMEC_L(5,:,0,2,2)&
                &-q_saf(:,0)*VMEC_L(5,:,1,2,2)
            T_loc(:,3,2) = 0.0_dp
            T_loc(:,3,3) = VMEC_L(5,:,1,2,2)
        end subroutine
    end subroutine
    
    subroutine test_calc_inv_met
        use metric_ops, only: calc_inv_met, &
            &T_VF, T_FV
        use eq_ops, only: calc_eq
        use eq_vars, only: n_par
        use VMEC_vars, only: n_r
        
        ! local variables
        integer :: id, kd                                                       ! counters
        real(dp) :: TxT(3,3)
        real(dp) :: trace(n_par,n_r)
        
        call writo('test calc_inv_met?')
        if(yes_no(.false.)) then
            call writo('Testing whether T_VF is equal to T_FV')
            call lvl_ud(1)
            
            call writo('calculating equilibrium')
            call lvl_ud(1)
            call calc_eq(0.2*pi)
            call lvl_ud(-1)
            
            call writo('testing zeroth order')
            call lvl_ud(1)
            do kd = 1,n_r
                do id = 1,n_par
                    TxT = matmul(T_VF(id,kd,:,:,0,0,0),T_FV(id,kd,:,:,0,0,0))
                    trace(id,kd) = TxT(1,1)+TxT(2,2)+TxT(3,3)
                end do
            end do
            call writo('trace of T_VF*T_FV = ')
            call print_ar_2(trace)
            call lvl_ud(-1)
            read(*,*)
            
            call writo('testing first order: r derivative')
            call lvl_ud(1)
            do kd = 1,n_r
                do id = 1,n_par
                    TxT = matmul(T_VF(id,kd,:,:,1,0,0),T_FV(id,kd,:,:,0,0,0)) &
                        &+ matmul(T_VF(id,kd,:,:,0,0,0),T_FV(id,kd,:,:,1,0,0))
                    trace(id,kd) = sum(TxT)
                end do
            end do
            call writo('norm of Dr(T_VF*T_FV) = ')
            call print_ar_2(trace)
            call lvl_ud(-1)
            read(*,*)
            
            call writo('testing first order: t derivative')
            call lvl_ud(1)
            do kd = 1,n_r
                do id = 1,n_par
                    TxT = matmul(T_VF(id,kd,:,:,0,1,0),T_FV(id,kd,:,:,0,0,0)) &
                        &+ matmul(T_VF(id,kd,:,:,0,0,0),T_FV(id,kd,:,:,0,1,0))
                    trace(id,kd) = sum(TxT)
                end do
            end do
            call writo('norm of Dt(T_VF*T_FV) = ')
            call print_ar_2(trace)
            call lvl_ud(-1)
            read(*,*)
            
            call writo('testing first order: z derivative')
            call lvl_ud(1)
            do kd = 1,n_r
                do id = 1,n_par
                    TxT = matmul(T_VF(id,kd,:,:,0,0,1),T_FV(id,kd,:,:,0,0,0)) &
                        &+ matmul(T_VF(id,kd,:,:,0,0,0),T_FV(id,kd,:,:,0,0,1))
                    trace(id,kd) = sum(TxT)
                end do
            end do
            call writo('norm of Dz(T_VF*T_FV) = ')
            call print_ar_2(trace)
            call lvl_ud(-1)
            read(*,*)
            
            call writo('testing second order: tt derivative')
            call lvl_ud(1)
            do kd = 1,n_r
                do id = 1,n_par
                    TxT = matmul(T_VF(id,kd,:,:,0,2,0),T_FV(id,kd,:,:,0,0,0)) &
                        &+2*matmul(T_VF(id,kd,:,:,0,1,0),T_FV(id,kd,:,:,0,1,0))&
                        &+matmul(T_VF(id,kd,:,:,0,0,0),T_FV(id,kd,:,:,0,2,0))
                    trace(id,kd) = sum(TxT)
                end do
            end do
            call writo('norm of Dtt(T_VF*T_FV) = ')
            call print_ar_2(trace)
            call lvl_ud(-1)
            read(*,*)
            
            call writo('testing third order: rtt derivative')
            call lvl_ud(1)
            do kd = 1,n_r
                do id = 1,n_par
                    TxT = matmul(T_VF(id,kd,:,:,1,2,0),T_FV(id,kd,:,:,0,0,0)) &
                        &+matmul(T_VF(id,kd,:,:,0,2,0),T_FV(id,kd,:,:,1,0,0))&
                        &+2*matmul(T_VF(id,kd,:,:,1,1,0),T_FV(id,kd,:,:,0,1,0))&
                        &+2*matmul(T_VF(id,kd,:,:,0,1,0),T_FV(id,kd,:,:,1,1,0))&
                        &+matmul(T_VF(id,kd,:,:,1,0,0),T_FV(id,kd,:,:,0,2,0))&
                        &+matmul(T_VF(id,kd,:,:,0,0,0),T_FV(id,kd,:,:,1,2,0))
                    trace(id,kd) = sum(TxT)
                end do
            end do
            call writo('norm of Drtt(T_VF*T_FV) = ')
            call print_ar_2(trace)
            call lvl_ud(-1)
            read(*,*)
            
            call lvl_ud(-1)
        end if
    end subroutine

    subroutine test_det
        use utilities, only: det
        
        ! local variables
        integer :: test_nr
        
        call writo('test det?')
        if(yes_no(.false.)) then
            call writo('Testing whether det is correctly calculated')
            call lvl_ud(1)
            do
                call writo('Select 1: 0D or 2: 2D')
                read(*,*) test_nr
                if (test_nr.gt.0 .and. test_nr.lt.3) then
                    exit
                else
                    call writo('You have to make a correct selection')
                end if
            end do
            call lvl_ud(1)
            
            select case (test_nr)
                case (1)
                    call test_0
                case (2)
                    call test_2
            end select
            
            call lvl_ud(-1)
            call lvl_ud(-1)
        end if
    contains
        subroutine test_0
            ! local variables
            real(dp) :: A(3,3)
            real(dp) :: detA
            
            ! define a 3x3 matrix
            A = transpose(reshape([1,2,3,0,3,5,7,8,0],[3,3]))
            call writo('calculating determinant of:')
            call print_ar_2(A)
            
            detA = 0.0_dp
            detA = det(A)
            call writo('det(A) = '//trim(r2str(detA)))
            call writo('compare with analyitcal result: -33')
        end subroutine
        
        subroutine test_2
            ! local variables
            real(dp), allocatable :: A(:,:,:,:)
            real(dp), allocatable :: detA(:,:), detA_ALT(:,:)
            integer :: id, jd, kd                                               ! counters
            real(dp), allocatable :: A_part(:,:)
            
            ! initialize
            allocate(A(10,5,4,4)); A = 0.0_dp
            allocate(detA(10,5)); detA = 0.0_dp
            allocate(detA_ALT(10,5)); detA_ALT = 0.0_dp
            allocate(A_part(3,3)); A_part = 0.0_dp
            
            !  define  4x4  matrix  over   2D  field  and  calculate  individual
            ! determinants
            do kd = 1,5
                do id = 1,10
                    A(id,kd,:,:) = transpose(reshape(&
                        &[(1.0_dp+kd*id**2/sqrt(jd*1.0_dp), jd=1,16)],[4,4]))
                    A_part = A(id,kd,:,:)
                    detA(id,kd) = det(A_part)
                end do
            end do
            
            ! calculate alternative all at once 
            detA_ALT = det(A)
            
            ! difference
            call writo('determinants calculated individually:')
            call print_ar_2(detA)
            call writo('determinants calculated global:')
            call print_ar_2(detA_ALT)
            call writo('difference')
            call print_ar_2(detA-detA_ALT)
        end subroutine
    end subroutine
    
    subroutine test_inv()
        use utilities, only: inv
        
        ! local variables
        real(dp), allocatable :: A(:,:,:,:), A_inv(:,:,:,:)
        integer :: id, jd, kd                                                   ! counters
        
        call writo('test inv?')
        if(yes_no(.false.)) then
            call writo('Testing whether inv works correctly')
            call lvl_ud(1)
            
            ! initialize
            allocate(A(2,5,4,4)); A = 0.0_dp
            allocate(A_inv(2,5,4,4)); A_inv = 0.0_dp
            
            ! define  4x4   matrix  over  2D  field   and  calculate  individual
            ! determinants
            do kd = 1,5
                do id = 1,2
                    A(id,kd,:,:) = transpose(reshape(&
                        &[(1.0_dp+kd*id**2/sqrt(jd*1.0_dp), jd=1,16)],[4,4]))
                end do
            end do
            
            ! calculate inverse
            A_inv = A
            A_inv = inv(A)
            
            ! multiply
            do kd = 1,5
                do id = 1,2
                    call writo('(id,kd) = ('//trim(i2str(id))//','&
                        &//trim(i2str(kd))//')')
                    call writo('matrix A = ')
                    call print_ar_2(A(id,kd,:,:))
                    call writo('inverse of A = ')
                    call print_ar_2(A_inv(id,kd,:,:))
                    call writo('product = ')
                    call print_ar_2(matmul(A(id,kd,:,:),A_inv(id,kd,:,:)))
                    call writo('')
                    read(*,*)
                end do
            end do
            call lvl_ud(-1)
        end if
    end subroutine
    
    subroutine test_metric_transf
        use eq_ops, only: calc_eq
        use eq_vars, only: n_par, flux_p, VMEC_R, VMEC_Z, VMEC_L, theta, zeta
        use VMEC_vars, only: n_r, B_V_s_H, B_V_c_H, mpol, ntor
        use fourier_ops, only: mesh_cs, f2r
        use metric_ops, only: jac_V, jac_F, T_FV, det_T_FV, g_F
        use utilities, only: VMEC_norm_deriv, det
        
        ! local variables
        real(dp) :: jac_ALT(1:n_par,1:n_r)
        integer :: id, kd                                                       ! counter
        real(dp) :: B(1:n_par,1:n_r)                                            ! magn. field
        real(dp) :: B_ALT(1:n_par,1:n_r)                                        ! magn. field
        real(dp) :: cs(0:mpol-1,-ntor:ntor,2)                                   ! (co)sines for all pol m and tor n
        
        call writo('test metric_transf?')
        if(yes_no(.false.)) then
            call lvl_ud(1)
            
            call writo('calculating equilibrium')
            call lvl_ud(1)
            call calc_eq(0.2*pi)
            call lvl_ud(-1)
            
            call writo('checking J_V')
            !   jac_F = (flux_p'/2pi * (1+L_t))^-1 * R * (R'Z_t - R' Z_t)
            call lvl_ud(1)
            do kd = 1,n_r
                jac_ALT(:,kd) = VMEC_R(:,kd,0,0,0)*(VMEC_R(:,kd,1,0,0)*&
                    &VMEC_Z(:,kd,0,1,0)-VMEC_R(:,kd,0,1,0)*VMEC_Z(:,kd,1,0,0))
            end do
            call writo('jac_V =')
            call print_ar_2(jac_V(:,:,0,0,0))
            call writo('jac_V_ALT =')
            call print_ar_2(jac_ALT)
            call writo('jac_V-jac_V_ALT =')
            call print_ar_2(jac_V(:,:,0,0,0)-jac_ALT)
            read(*,*)
            call lvl_ud(-1)
            
            call writo('checking J_F')
            !   jac_F = (flux_p'/2pi * (1+L_t))^-1 * R * (R'Z_t - R' Z_t)
            call lvl_ud(1)
            do kd = 1,n_r
                jac_ALT(:,kd) = VMEC_R(:,kd,0,0,0)*(VMEC_R(:,kd,1,0,0)*&
                    &VMEC_Z(:,kd,0,1,0)-VMEC_R(:,kd,0,1,0)*VMEC_Z(:,kd,1,0,0))&
                    &/(flux_p(kd,1)/(2*pi)*(1+VMEC_L(:,kd,0,1,0)))
            end do
            call writo('jac_F =')
            call print_ar_2(jac_F(:,:,0,0,0))
            call writo('jac_F_ALT =')
            call print_ar_2(jac_ALT)
            call writo('jac_F-jac_F_ALT =')
            call print_ar_2(jac_F(:,:,0,0,0)-jac_ALT)
            read(*,*)
            call lvl_ud(-1)
            
            call writo('checking r derivatives')
            call write_out(1,n_r,jac_F(5,:,0,0,0),'J_F(par=5)')
            call VMEC_norm_deriv(jac_F(5,:,0,0,0),jac_ALT(5,:),n_r-1._dp,1,1)
            call write_out(1,n_r,jac_F(5,:,1,0,0),'Dr J_F(par=5)')
            call write_out(1,n_r,jac_ALT(5,:)-jac_F(5,:,1,0,0),&
                &'diff with num Dr J_F(par=5)')
            call writo('checking rr derivatives')
            call VMEC_norm_deriv(jac_F(5,:,0,0,0),jac_ALT(5,:),n_r-1._dp,2,1)
            call write_out(1,n_r,jac_F(5,:,3,0,0),'Drr J_F(par=5)')
            call write_out(1,n_r,jac_ALT(5,:)-jac_F(5,:,2,0,0),&
                &'diff with num Drr J_F(par=5)')
            
            call writo('checking jacobians as determinants')
            call lvl_ud(1)
            jac_ALT = det(T_FV(:,:,:,:,0,0,0))
            call print_ar_2(jac_ALT-det_T_FV(:,:,0,0,0))
            read(*,*)
            call lvl_ud(-1)
            
            call writo('checking magnetic field')
            call lvl_ud(1)
            do kd = 1,n_r
                do id = 1,n_par
                    cs = mesh_cs(mpol,ntor,theta(id,kd),zeta(id,kd))
                    B(id,kd) = &
                        &f2r(B_V_c_H(:,:,kd),B_V_s_H(:,:,kd),cs,mpol,ntor,[0,0])
                    B_ALT(id,kd) = sqrt(g_F(id,kd,2,2,0,0,0))/jac_F(id,kd,0,0,0)
                end do
            end do
            call writo('B = ')
            call print_ar_2(B)
            call writo('alternative B = ')
            call print_ar_2(B_ALT)
            call write_out(1,n_r,B(5,:),'B(par=5)')
            call write_out(1,n_r-1,B_ALT(5,2:n_r),'B_ALT(par=5)')
            call lvl_ud(-1)
            
            call lvl_ud(-1)
        end if
    end subroutine

    subroutine test_B
        !use eq_ops, only: calc_eq
        !use eq_vars, only: n_par, flux_p, flux_p_H, q_saf, q_saf_H, lam_H
        !use metric_ops, only: jac_V, jac_V_H, g_V, g_V_H, jac_F, jac_F_H
        !use VMEC_vars, only: n_r
        !use utilities, only: h2f
        !use B_vars, only: B_V_sub, B_V_sub_H, B_V, B_V_H
        
        !! local variables
        !real(dp) :: B_V_sub_alt(n_par,n_r,3)
        !real(dp) :: B_V_sub_alt2(n_par,n_r,3)
        !real(dp) :: B_V_sub_H_alt(n_par,n_r,3)
        !real(dp) :: B_V_sub_H_alt2(n_par,n_r,3)
        !real(dp) :: B_V_alt(n_par,n_r)
        !real(dp) :: B_V_H_alt(n_par,n_r)
        !real(dp) :: diff(n_r)
        !integer :: id, jd, kd
        !real(dp) :: lam(n_par,n_r,4)
        
        !call writo('test magnetic fields?')
        !if(yes_no(.false.)) then
            !call lvl_ud(1)
            
            !! calculate equilibrium with some value for alpha
            !call writo('Calculate equilibrium')
            !call lvl_ud(1)
            !call calc_eq(pi/2)
            !call lvl_ud(-1)
            
            !call writo('Testing whether the magnetic fields from VMEC are &
                !&in agreement with B_i = nabla alpha x nabla psi')
            !call lvl_ud(1)
            
            !! calculate lam from lam_H
            !do id = 1,n_par
                !do jd = 1,4
                    !lam(id,:,jd) = h2f(lam_H(id,:,jd))
                !end do
            !end do
            
            !! calculate B_i alternatively
            !B_V_sub_alt = 0.0_dp
            !B_V_sub_H_alt = 0.0_dp
            !perp: do kd = 1,n_r                                                 ! avoid problems at magn. axis.
                !par: do id = 1,n_par
                    !comp: do jd = 1,3
                        !B_V_sub_alt(id,kd,jd) = &
                            !&flux_p(kd,2)/(2*pi*jac_V(id,kd)) * &
                            !&(q_saf(kd,1)*(1+lam(id,kd,3)) * g_V(3,jd,id,kd) + &
                            !&(1-q_saf(kd,1)*lam(id,kd,4)) * g_V(2,jd,id,kd))
                        !B_V_sub_alt2(id,kd,jd) = &
                            !&1/jac_F(id,kd) * &
                            !&(q_saf(kd,1) * g_V(3,jd,id,kd) + &
                            !&(1-q_saf(kd,1)*lam(id,kd,4))/(1+lam(id,kd,3))* &
                            !&g_V(2,jd,id,kd))
                        !B_V_sub_H_alt(id,kd,jd) = &
                            !&flux_p_H(kd,2)/(2*pi*jac_V_H(id,kd)) * &
                            !&(q_saf_H(kd,1)*(1+lam_H(id,kd,3)) &
                            !&* g_V_H(3,jd,id,kd) + (1-q_saf_H(kd,1)*&
                            !&lam_H(id,kd,4)) * g_V_H(2,jd,id,kd))
                        !B_V_sub_H_alt2(id,kd,jd) = &
                            !&1/jac_F_H(id,kd) * &
                            !&(q_saf_H(kd,1) * g_V_H(3,jd,id,kd) + &
                            !&(1-q_saf_H(kd,1)*lam_H(id,kd,4))/&
                            !&(1+lam_H(id,kd,3))* g_V_H(2,jd,id,kd))
                    !end do comp
                !end do par
            !end do perp
            
            !! calculate B alternatively
            !B_V_alt = 0.0_dp
            !B_V_H_alt = 0.0_dp
            !perp2: do kd = 1,n_r
                !par2: do id = 1,n_par
                    !B_V_alt(id,kd) = sqrt(1/(jac_F(id,kd)) * &
                        !&((1-q_saf(kd,1)*lam(id,kd,4))/(1+lam(id,kd,3)) &
                        !&* B_V_sub(id,kd,1,2) + q_saf(kd,1)*B_V_sub(id,kd,1,3)))
                    !B_V_H_alt(id,kd) = sqrt(1/(jac_F_H(id,kd)) * &
                        !&((1-q_saf_H(kd,1)*lam_H(id,kd,4))/(1+lam_H(id,kd,3)) &
                        !&* B_V_sub_H(id,kd,1,2) &
                        !&+ q_saf_H(kd,1)*B_V_sub_H(id,kd,1,3)))
                !end do par2
            !end do perp2
            
            !call writo('FULL MESH')
            !call lvl_ud(1)
            !diff = 0.0_dp
            !call writo('Paused... press enter')
            !read(*,*)
            
            !! FM magnitude
            !do kd = 1, n_r
                !diff(kd) = log10(min(max(maxval(abs(B_V(:,kd)-B_V_alt(:,kd)))&
                    !&,1.0E-20_dp),1E10_dp))
            !end do
            !call write_out(1,n_r,diff,'test whether sqrt(B_theta/J)_F = B &
                !&[log]')
            
            !! FM r component
            !do kd = 1, n_r
                !diff(kd) = log10(min(max(maxval(abs(B_V_sub(:,kd,1,1)-&
                    !&B_V_sub_alt(:,kd,1))),1.0E-20_dp),1E10_dp))
            !end do
            !call write_out(1,n_r,diff,'max diff. between standard and &
                !&alternative treatment, r component [log]')
            !do kd = 1, n_r
                !diff(kd) = log10(min(max(maxval(abs(B_V_sub(:,kd,1,1)-&
                    !&B_V_sub_alt2(:,kd,1))),1.0E-20_dp),1E10_dp))
            !end do
            !call write_out(1,n_r,diff,'max diff. between standard and &
                !&alternative treatment 2, r component [log]')
            
            !! FM theta component
            !do kd = 1, n_r
                !diff(kd) = log10(min(max(maxval(abs(B_V_sub(:,kd,1,2)-&
                    !&B_V_sub_alt(:,kd,2))),1.0E-20_dp),1E10_dp))
            !end do
            !call write_out(1,n_r,diff,'max diff. between standard and &
                !&alternative treatment, theta component [log]')
            !do kd = 1, n_r
                !diff(kd) = log10(min(max(maxval(abs(B_V_sub(:,kd,1,2)-&
                    !&B_V_sub_alt2(:,kd,2))),1.0E-20_dp),1E10_dp))
            !end do
            !call write_out(1,n_r,diff,'max diff. between standard and &
                !&alternative treatment 2, theta component [log]')
                
            !! FM zeta component
            !do kd = 1, n_r
                !diff(kd) = log10(min(max(maxval(abs(B_V_sub(:,kd,1,3)-&
                    !&B_V_sub_alt(:,kd,3))),1.0E-20_dp),1E10_dp))
            !end do
            !call write_out(1,n_r,diff,'max diff. between standard and &
                !&alternative treatment, zeta component [log]')
            !do kd = 1, n_r
                !diff(kd) = log10(min(max(maxval(abs(B_V_sub(:,kd,1,3)-&
                    !&B_V_sub_alt2(:,kd,3))),1.0E-20_dp),1E10_dp))
            !end do
            !call write_out(1,n_r,diff,'max diff. between standard and &
                !&alternative treatment 2, zeta component [log]')
            
            !call lvl_ud(-1)
            
            !call writo('HALF MESH')
            !call lvl_ud(1)
            !call writo('Paused... press enter')
            !read(*,*)
            
            !! HM magnitude
            !do kd = 1, n_r
                !diff(kd) = log10(min(max(maxval(abs(B_V_H(:,kd)-&
                    !&B_V_H_alt(:,kd))),1.0E-20_dp),1E10_dp))
            !end do
            !call write_out(1,n_r-1,diff(2:n_r),'test whether &
                !&sqrt(B_theta/J)_F = B [log]')
            
            !! HM r component
            !do kd = 1, n_r
                !diff(kd) = log10(max(maxval(abs(B_V_sub_H(:,kd,1,1)-&
                    !&B_V_sub_H_alt(:,kd,1))),1.0E-20_dp))
            !end do
            !call write_out(1,n_r-1,diff(2:n_r),'max diff. between standard and &
                !&alternative treatment, r component [log]')
            !do kd = 1, n_r
                !diff(kd) = log10(max(maxval(abs(B_V_sub_H(:,kd,1,1)-&
                    !&B_V_sub_H_alt2(:,kd,1))),1.0E-20_dp))
            !end do
            !call write_out(1,n_r-1,diff(2:n_r),'max diff. between standard and &
                !&alternative treatment 2, r component [log]')
            
            !! HM theta component
            !do kd = 1, n_r
                !diff(kd) = log10(max(maxval(abs(B_V_sub_H(:,kd,1,2)-&
                    !&B_V_sub_H_alt(:,kd,2))),1.0E-20_dp))
            !end do
            !call write_out(1,n_r-1,diff(2:n_r),'max diff. between standard and &
                !&alternative treatment, theta component [log]')
            !do kd = 1, n_r
                !diff(kd) = log10(max(maxval(abs(B_V_sub_H(:,kd,1,2)-&
                    !&B_V_sub_H_alt2(:,kd,2))),1.0E-20_dp))
            !end do
            !call write_out(1,n_r-1,diff(2:n_r),'max diff. between standard and &
                !&alternative treatment 2, theta component [log]')
            
            !! HM zeta component
            !do kd = 1, n_r
                !diff(kd) = log10(max(maxval(abs(B_V_sub_H(:,kd,1,3)-&
                    !&B_V_sub_H_alt(:,kd,3))),1.0E-20_dp))
            !end do
            !call write_out(1,n_r-1,diff(2:n_r),'max diff. between standard and &
                !&alternative treatment, zeta component [log]')
            !do kd = 1, n_r
                !diff(kd) = log10(max(maxval(abs(B_V_sub_H(:,kd,1,3)-&
                    !&B_V_sub_H_alt2(:,kd,3))),1.0E-20_dp))
            !end do
            !call write_out(1,n_r-1,diff(2:n_r),'max diff. between standard and &
                !&alternative treatment 2, zeta component [log]')
            
            !call lvl_ud(-1)
            
            !call lvl_ud(-1)
            !call lvl_ud(-1)
        !end if
    end subroutine

    subroutine test_ang_B
        use VMEC_vars, only: n_r, mpol, ntor
        use eq_vars, only: eqd_mesh, calc_mesh, &
            &n_par, theta, zeta
        use output_ops, only: format_out
        use num_vars, only: theta_var_along_B
        
        ! local variables (not to be used in child functions)
        real(dp) :: alpha
        integer :: id
        real(dp), allocatable :: plot_ang(:,:)
        real(dp), allocatable :: plot_dep_var(:)
        real(dp), allocatable :: plot_var(:,:)
        integer :: n_plot
        character(len=max_str_ln) :: par_ang, dep_ang                           ! parallel angle and dependent angle, for output message
        integer :: format_out_old
        real(dp) :: grid_min, grid_max
        integer :: n_par_old
        
        ! local variables (also used in child functions)
        integer :: kd
        
        call writo('test ang_B?')
        if(yes_no(.false.)) then
            call writo('Plot zeta(theta)')
            
            call lvl_ud(1)
            
            ! CALCULATE THETA(ZETA)
            call writo('Calculating theta(zeta)')
            alpha = 1.2_dp*pi
            
            ! decrease the number of parallel points for less plots
            n_par_old = n_par
            !n_par = 4
            
            ! calculate mesh points (theta, zeta) that follow the magnetic field
            ! line
            call calc_mesh(alpha)
            
            allocate(plot_ang(2,n_par))
            do kd = 1,n_r
                ! plot (theta,zeta)
                plot_ang(1,:) = theta(:,kd)
                plot_ang(2,:) = zeta(:,kd)
                
                ! plot it on the screen using format_out = 3
                format_out_old = format_out
                format_out = 3
                call write_out(2, n_par, plot_ang, &
                        &'zeta(theta) at r = '//trim(i2str(kd))&
                        &//'/'//trim(i2str(n_r)),comment=&
                        &trim(r2strt(minval(plot_ang(1,:))))//' < theta < '&
                        &//trim(r2strt(maxval(plot_ang(1,:))))//' and '//&
                        &trim(r2strt(minval(plot_ang(2,:))))//' < zeta < '&
                        &//trim(r2strt(maxval(plot_ang(2,:)))))
                format_out = format_out_old
                
                call writo('Paused... plot next?')
                if(.not.yes_no(.true.)) then
                    exit
                end if
            end do
            
            call lvl_ud(-1)
            
            ! CALCULATE F FOR A RANGE OF PARALLEL VALUES
            call writo('Calculating f for a range of parallel values')
            
            call lvl_ud(1)
            n_plot = 100
            
            ! initialize
            allocate(plot_dep_var(n_plot)); plot_dep_var = 0.0_dp
            allocate(plot_var(2,n_plot)); plot_var = 0.0_dp
            
            ! set correct plot messages
            if (theta_var_along_B) then                                         ! looking for zeta
                par_ang = 'theta'; dep_ang = 'zeta'
            else                                                                ! looking for theta
                par_ang = 'zeta'; dep_ang = 'theta'
            end if
                
            perp: do kd = 1, n_r
                ! determine grid considering the solutions on current flux surface
                if (theta_var_along_B) then                                         ! looking for zeta
                    grid_min = minval(zeta(:,kd)) - &
                        &0.1*(maxval(zeta(:,kd))-minval(zeta(:,kd)))
                    grid_max = maxval(zeta(:,kd)) + &
                        &0.1*(maxval(zeta(:,kd))-minval(zeta(:,kd)))
                else                                                                ! looking for theta
                    grid_min = minval(theta(:,kd)) - &
                        &0.1*(maxval(theta(:,kd))-minval(theta(:,kd)))
                    grid_max = maxval(theta(:,kd)) + &
                        &0.1*(maxval(theta(:,kd))-minval(theta(:,kd)))
                end if
                plot_dep_var = eqd_mesh(n_plot, grid_min, grid_max)
                par: do id = 1, n_par
                    if (theta_var_along_B) then                                 ! looking for zeta
                        plot_var(1,:) = plot_dep_var                            ! eq. mesh of zeta
                        plot_var(2,:) = find_f_plot(n_plot,plot_dep_var, &
                            &theta(id,kd),alpha)                              ! corresponding f(zeta)
                    else                                                        ! looking for theta
                        plot_var(1,:) = plot_dep_var                            ! eq. mesh of theta
                        plot_var(2,:) = find_f_plot(n_plot,plot_dep_var, &
                            &zeta(id,kd),alpha)                               ! corresponding f(theta)
                    end if
                    
                    ! plot it on the screen using format_out = 3
                    format_out_old = format_out
                    format_out = 3
                    call write_out(2,n_plot,plot_var, 'f('//trim(dep_ang)&
                        &//') = zeta -q(theta + lambda) - alpha_0 at r = '//&
                        &trim(i2str(kd))//'/'//trim(i2str(n_r)), comment=&
                        &'= 0 at (theta, zeta) = ('//&
                        &trim(r2strt(theta(id,kd)))&
                        &//', '//trim(r2strt(zeta(id,kd)))//')')
                    format_out = format_out_old
                    
                    call writo('Paused... plot next?')
                    if(.not.yes_no(.true.)) then
                        exit perp
                    end if
                end do par
            end do perp
            
            n_par = n_par_old
            
            call lvl_ud(-1)
        end if
    contains
        ! calculates the function f = zeta - q (theta + lambda) - alpha_0
        ! makes use of kd from parent to indicate the flux surface
        function find_f_plot(n_ang,dep_var,par_var,alpha)
            use fourier_ops, only: mesh_cs, f2r
            use VMEC_vars, only: iotaf
            
            ! input / output
            integer :: n_ang
            real(dp) :: find_f_plot(n_ang)
            real(dp) :: dep_var(n_ang)
            real(dp) :: par_var
            real(dp) :: alpha
            
            ! local variables
            integer :: jd
            real(dp), allocatable :: cs(:,:,:)                                  ! (co)sines for all pol m and tor n
            real(dp) :: lam
            real(dp) :: l_c_F(0:mpol-1,-ntor:ntor,1:n_r)                        ! FM version of HM l_c
            real(dp) :: l_s_F(0:mpol-1,-ntor:ntor,1:n_r)                        ! FM version of HM l_s
            
            allocate(cs(0:mpol-1,-ntor:ntor,2))
            do jd = 1, n_ang
                if (theta_var_along_B) then                                     ! looking for zeta (dep. var)
                    ! cosines and sines
                    cs = mesh_cs(mpol,ntor,par_var,dep_var(jd))
                    ! lambda
                    lam = f2r(l_c_F(:,:,kd),l_s_F(:,:,kd),cs,mpol,ntor)
                    find_f_plot(jd) = dep_var(jd)-(par_var+lam)/iotaf(kd)-alpha
                else                                                            ! looking for theta (dep. var)
                    ! cosines and sines
                    cs = mesh_cs(mpol,ntor,dep_var(jd),par_var)
                    ! lambda
                    lam = f2r(l_c_F(:,:,kd),l_s_F(:,:,kd),cs,mpol,ntor)
                    find_f_plot(jd) = par_var-(dep_var(jd)+lam)/iotaf(kd)-alpha
                end if
            end do
        end function find_f_plot
    end subroutine
    
    subroutine test_arr_mult
        use utilities, only: arr_mult
        use eq_vars, only: eqd_mesh, calc_RZL, init_eq, calc_flux_q, &
            &VMEC_R, vmec_Z, n_par, theta, zeta, q_saf, pres
        use VMEC_vars, only: n_r
        !use num_vars, only: max_deriv
        
        ! local variables
        integer :: case_nr
        
        call writo('test arr_mult?')
        if(yes_no(.false.)) then
            call writo('Testing whether arr_mult correctly multiplies two &
                &functions, taking into account derivatives')
            call lvl_ud(1)
            do 
                case_nr = 1
                call writo('Choose which case to test')
                call writo('    1. arr_1 in 3 coords and arr_2 in 3 coords [def]')
                call writo('    2. arr_1 in 3 coords and arr_2 in 1 coords')
                call writo('    2. arr_1 in 1 coords and arr_2 in 1 coords')
                read(*,*) case_nr
                if (case_nr.eq.1 .or. case_nr.eq.2 .or. case_nr.eq.3) then
                    exit
                else
                    call writo('You have to choose one of the three &
                        &possibilities...')
                end if
            end do
            
            call lvl_ud(1)
            select case (case_nr)
                case(1)
                    call case_3_3
                case(2)
                    call case_3_1
                case(3)
                    call case_1_1
                case default
                    call writo('ERROR: no case number associated with '&
                        &//trim(i2str(case_nr)))
                    call writo('How did you get here???')
                    stop
            end select
            
            call lvl_ud(-1)
            call lvl_ud(-1)
        end if
    contains
        subroutine case_3_3
            ! local variables
            real(dp) :: RZ(n_par,n_r), RZ_num(n_r)
            integer :: id, jd, kd
            
            ! calculate equilibrium quantities
            call writo('calculate RZL')
            call lvl_ud(1)
            if (allocated(zeta)) deallocate(zeta)
            allocate(zeta(n_par,n_r)); zeta = 0.0_dp
            if (allocated(theta)) deallocate(theta)
            allocate(theta(n_par,n_r)); theta = 0.0_dp
            
            call init_eq
            
            zeta = 0.4*pi/2
            do jd = 1,n_r
                theta(:,jd) = eqd_mesh(n_par, 0.0_dp*pi, 3.0_dp*pi)
            end do
            
            do kd = 0, 3
                do jd = 0, 3
                    do id = 0, 3
                        call calc_RZL([id,jd,kd])
                    end do
                end do
            end do
            call lvl_ud(-1)
            
            ! multiply
            call writo('multiply R and Z')
            call lvl_ud(1)
            RZ = 0.0_dp
            call arr_mult(VMEC_R,VMEC_Z,RZ,[0,0,0])
            RZ_num = VMEC_R(5,:,0,0,0)*VMEC_Z(5,:,0,0,0)
            call write_out(1,n_r,RZ(5,:),'RZ at par = 5')
            call write_out(1,n_r,RZ_num,'check on RZ at par = 5',comment=&
                &'max diff. = '//trim(r2str(maxval(abs(RZ(5,:)-RZ_num)))))
            call lvl_ud(-1)
            
            ! derive multiplied values
            call writo('derive RZ')
            call lvl_ud(1)
            ! Dr
            RZ = 0.0_dp
            call arr_mult(VMEC_R,VMEC_Z,RZ,[1,0,0])
            RZ_num = VMEC_R(5,:,1,0,0)*VMEC_Z(5,:,0,0,0) + &
                &VMEC_R(5,:,0,0,0)*VMEC_Z(5,:,1,0,0)
            call write_out(1,n_r,RZ(5,:),'DrRZ at par = 5')
            call write_out(1,n_r,RZ_num,'check on DrRZ at par = 5',comment=&
                &'max diff. = '//trim(r2str(maxval(abs(RZ(5,:)-RZ_num)))))
            ! Dtheta
            RZ = 0.0_dp
            call arr_mult(VMEC_R,VMEC_Z,RZ,[0,1,0])
            RZ_num = VMEC_R(5,:,0,1,0)*VMEC_Z(5,:,0,0,0) + &
                &VMEC_R(5,:,0,0,0)*VMEC_Z(5,:,0,1,0)
            call write_out(1,n_r,RZ(5,:),'DtRZ at par = 5')
            call write_out(1,n_r,RZ_num,'check on DtRZ at par = 5',comment=&
                &'max diff. = '//trim(r2str(maxval(abs(RZ(5,:)-RZ_num)))))
            ! Dzeta
            RZ = 0.0_dp
            call arr_mult(VMEC_R,VMEC_Z,RZ,[0,0,1])
            RZ_num = VMEC_R(5,:,0,0,1)*VMEC_Z(5,:,0,0,0) + &
                &VMEC_R(5,:,0,0,0)*VMEC_Z(5,:,0,0,1)
            call write_out(1,n_r,RZ(5,:),'DzRZ at par = 5')
            call write_out(1,n_r,RZ_num,'check on DzRZ at par = 5',comment=&
                &'max diff. = '//trim(r2str(maxval(abs(RZ(5,:)-RZ_num)))))
            call lvl_ud(-1)
            
            ! double derive multiplied values
            call writo('double derive RZ')
            call lvl_ud(1)
            ! Drr
            RZ = 0.0_dp
            call arr_mult(VMEC_R,VMEC_Z,RZ,[2,0,0])
            call write_out(1,n_r,RZ(5,:),'DrrRZ at par = 5')
            call write_out(1,n_r,VMEC_R(5,:,2,0,0)*VMEC_Z(5,:,0,0,0) + &
                &2.0_dp * VMEC_R(5,:,1,0,0)*VMEC_Z(5,:,1,0,0) + &
                &VMEC_R(5,:,0,0,0)*VMEC_Z(5,:,2,0,0),&
                &'check on DrrRZ at par = 5')
            ! Drtheta
            RZ = 0.0_dp
            call arr_mult(VMEC_R,VMEC_Z,RZ,[1,1,0])
            call write_out(1,n_r,RZ(5,:),'DrrRZ at par = 5')
            call write_out(1,n_r,VMEC_R(5,:,1,1,0)*VMEC_Z(5,:,0,0,0) + &
                &VMEC_R(5,:,1,0,0)*VMEC_Z(5,:,0,1,0) + &
                &VMEC_R(5,:,0,1,0)*VMEC_Z(5,:,1,0,0) + &
                &VMEC_R(5,:,0,0,0)*VMEC_Z(5,:,1,1,0),&
                &'check on DrtRZ at par = 5')
            ! Drzeta
            RZ = 0.0_dp
            call arr_mult(VMEC_R,VMEC_Z,RZ,[1,0,1])
            call write_out(1,n_r,RZ(5,:),'DrzRZ at par = 5')
            call write_out(1,n_r,VMEC_R(5,:,1,0,1)*VMEC_Z(5,:,0,0,0) + &
                &VMEC_R(5,:,1,0,0)*VMEC_Z(5,:,0,0,1) + &
                &VMEC_R(5,:,0,0,1)*VMEC_Z(5,:,1,0,0) + &
                &VMEC_R(5,:,0,0,0)*VMEC_Z(5,:,1,0,1),&
                &'check on DrzRZ at par = 5')
            ! Dthetatheta
            RZ = 0.0_dp
            call arr_mult(VMEC_R,VMEC_Z,RZ,[0,2,0])
            call write_out(1,n_r,RZ(5,:),'DttRZ at par = 5')
            call write_out(1,n_r,VMEC_R(5,:,0,2,0)*VMEC_Z(5,:,0,0,0) + &
                &2.0_dp * VMEC_R(5,:,0,1,0)*VMEC_Z(5,:,0,1,0) + &
                &VMEC_R(5,:,0,0,0)*VMEC_Z(5,:,0,2,0),&
                &'check on DttRZ at par = 5')
            ! Dtzeta
            RZ = 0.0_dp
            call arr_mult(VMEC_R,VMEC_Z,RZ,[0,1,1])
            call write_out(1,n_r,RZ(5,:),'DtzRZ at par = 5')
            call write_out(1,n_r,VMEC_R(5,:,0,1,1)*VMEC_Z(5,:,0,0,0) + &
                &VMEC_R(5,:,0,1,0)*VMEC_Z(5,:,0,0,1) + &
                &VMEC_R(5,:,0,0,1)*VMEC_Z(5,:,0,1,0) + &
                &VMEC_R(5,:,0,0,0)*VMEC_Z(5,:,0,1,1),&
                &'check on DtzRZ at par = 5')
            ! Dzetazeta
            RZ = 0.0_dp
            call arr_mult(VMEC_R,VMEC_Z,RZ,[0,0,2])
            call write_out(1,n_r,RZ(5,:),'DzzRZ at par = 5')
            call write_out(1,n_r,VMEC_R(5,:,0,0,2)*VMEC_Z(5,:,0,0,0) + &
                &2.0_dp * VMEC_R(5,:,0,0,1)*VMEC_Z(5,:,0,0,1) + &
                &VMEC_R(5,:,0,0,0)*VMEC_Z(5,:,0,0,2),&
                &'check on DzzRZ at par = 5')
            call lvl_ud(-1)
            
            ! higher order derive multiplied values
            call writo('higher order derive RZ')
            call lvl_ud(1)
            RZ = 0.0_dp
            call arr_mult(VMEC_R,VMEC_Z,RZ,[1,2,3])
            call write_out(1,n_r,RZ(5,:),'DrttzzzRZ at par = 5')
            call write_out(1,n_r,VMEC_R(5,:,1,2,3)*VMEC_Z(5,:,0,0,0) + &
                &3.0_dp * VMEC_R(5,:,1,2,2)*VMEC_Z(5,:,0,0,1) + &
                &3.0_dp * VMEC_R(5,:,1,2,1)*VMEC_Z(5,:,0,0,2) + &
                &1.0_dp * VMEC_R(5,:,1,2,0)*VMEC_Z(5,:,0,0,3) + &
                &2.0_dp * VMEC_R(5,:,1,1,3)*VMEC_Z(5,:,0,1,0) + &
                &6.0_dp * VMEC_R(5,:,1,1,2)*VMEC_Z(5,:,0,1,1) + &
                &6.0_dp * VMEC_R(5,:,1,1,1)*VMEC_Z(5,:,0,1,2) + &
                &2.0_dp * VMEC_R(5,:,1,1,0)*VMEC_Z(5,:,0,1,3) + &
                &1.0_dp * VMEC_R(5,:,1,0,3)*VMEC_Z(5,:,0,2,0) + &
                &3.0_dp * VMEC_R(5,:,1,0,2)*VMEC_Z(5,:,0,2,1) + &
                &3.0_dp * VMEC_R(5,:,1,0,1)*VMEC_Z(5,:,0,2,2) + &
                &1.0_dp * VMEC_R(5,:,1,0,0)*VMEC_Z(5,:,0,2,3) + &
                &1.0_dp * VMEC_R(5,:,0,2,3)*VMEC_Z(5,:,1,0,0) + &
                &3.0_dp * VMEC_R(5,:,0,2,2)*VMEC_Z(5,:,1,0,1) + &
                &3.0_dp * VMEC_R(5,:,0,2,1)*VMEC_Z(5,:,1,0,2) + &
                &1.0_dp * VMEC_R(5,:,0,2,0)*VMEC_Z(5,:,1,0,3) + &
                &2.0_dp * VMEC_R(5,:,0,1,3)*VMEC_Z(5,:,1,1,0) + &
                &6.0_dp * VMEC_R(5,:,0,1,2)*VMEC_Z(5,:,1,1,1) + &
                &6.0_dp * VMEC_R(5,:,0,1,1)*VMEC_Z(5,:,1,1,2) + &
                &2.0_dp * VMEC_R(5,:,0,1,0)*VMEC_Z(5,:,1,1,3) + &
                &1.0_dp * VMEC_R(5,:,0,0,3)*VMEC_Z(5,:,1,2,0) + &
                &3.0_dp * VMEC_R(5,:,0,0,2)*VMEC_Z(5,:,1,2,1) + &
                &3.0_dp * VMEC_R(5,:,0,0,1)*VMEC_Z(5,:,1,2,2) + &
                &1.0_dp * VMEC_R(5,:,0,0,0)*VMEC_Z(5,:,1,2,3) ,&
                &'check on DrttzzzRZ at par = 5')
            call lvl_ud(-1)
        end subroutine
        
        subroutine case_3_1
            ! local variables
            real(dp) :: Rq(n_par,n_r), Rq_num(n_r)
            integer :: id, jd, kd
            
            ! calculate equilibrium quantities
            call writo('calculate RZL')
            call lvl_ud(1)
            if (allocated(zeta)) deallocate(zeta)
            allocate(zeta(n_par,n_r)); zeta = 0.0_dp
            if (allocated(theta)) deallocate(theta)
            allocate(theta(n_par,n_r)); theta = 0.0_dp
            
            call init_eq
            
            zeta = 0.4*pi/2
            do jd = 1,n_r
                theta(:,jd) = eqd_mesh(n_par, 0.0_dp*pi, 3.0_dp*pi)
            end do
            
            do kd = 0, 3
                do jd = 0, 3
                    do id = 0, 3
                        call calc_RZL([id,jd,kd])
                    end do
                end do
            end do
            call lvl_ud(-1)
            
            call writo('calculate q_saf')
            call lvl_ud(1)
            call calc_flux_q
            call lvl_ud(-1)
            
            ! multiply
            call writo('multiply R and q_saf')
            call lvl_ud(1)
            Rq = 0.0_dp
            call arr_mult(VMEC_R,q_saf,Rq,[0,0,0])
            Rq_num = VMEC_R(5,:,0,0,0)*q_saf(:,0)
            call write_out(1,n_r,Rq(5,:),'Rq at par = 5')
            call write_out(1,n_r,Rq_num,'check on Rq at par = 5',comment=&
                &'max diff. = '//trim(r2str(maxval(abs(Rq(5,:)-Rq_num)))))
            call lvl_ud(-1)
            
            ! derive multiplied values
            call writo('derive Rq')
            call lvl_ud(1)
            ! Dr
            Rq = 0.0_dp
            call arr_mult(VMEC_R,q_saf,Rq,[1,0,0])
            Rq_num = VMEC_R(5,:,1,0,0)*q_saf(:,0) + &
                &VMEC_R(5,:,0,0,0)*q_saf(:,1)
            call write_out(1,n_r,Rq(5,:),'DrRq at par = 5')
            call write_out(1,n_r,Rq_num,'check on DrRq at par = 5',comment=&
                &'max diff. = '//trim(r2str(maxval(abs(Rq(5,:)-Rq_num)))))
            ! Dtheta
            Rq = 0.0_dp
            call arr_mult(VMEC_R,q_saf,Rq,[0,1,0])
            Rq_num = VMEC_R(5,:,0,1,0)*q_saf(:,0) + &
                &VMEC_R(5,:,0,0,0)*0.0_dp
            call write_out(1,n_r,Rq(5,:),'DtRq at par = 5')
            call write_out(1,n_r,Rq_num,'check on DtRq at par = 5',comment=&
                &'max diff. = '//trim(r2str(maxval(abs(Rq(5,:)-Rq_num)))))
            ! Dzeta
            Rq = 0.0_dp
            call arr_mult(VMEC_R,q_saf,Rq,[0,0,1])
            Rq_num = VMEC_R(5,:,0,0,1)*q_saf(:,0) + &
                &VMEC_R(5,:,0,0,0)*0.0_dp
            call write_out(1,n_r,Rq(5,:),'DzRq at par = 5')
            call write_out(1,n_r,Rq_num,'check on DzRq at par = 5',comment=&
                &'max diff. = '//trim(r2str(maxval(abs(Rq(5,:)-Rq_num)))))
            call lvl_ud(-1)
            
            ! double derive multiplied values
            call writo('double derive Rq')
            call lvl_ud(1)
            ! Drr
            Rq = 0.0_dp
            call arr_mult(VMEC_R,q_saf,Rq,[2,0,0])
            call write_out(1,n_r,Rq(5,:),'DrrRq at par = 5')
            call write_out(1,n_r,VMEC_R(5,:,2,0,0)*q_saf(:,0) + &
                &2.0_dp * VMEC_R(5,:,1,0,0)*q_saf(:,1) + &
                &VMEC_R(5,:,0,0,0)*q_saf(:,2),&
                &'check on DrrRq at par = 5')
            ! Drtheta
            Rq = 0.0_dp
            call arr_mult(VMEC_R,q_saf,Rq,[1,1,0])
            call write_out(1,n_r,Rq(5,:),'DrrRq at par = 5')
            call write_out(1,n_r,VMEC_R(5,:,1,1,0)*q_saf(:,0) + &
                &VMEC_R(5,:,1,0,0)*0.0_dp + &
                &VMEC_R(5,:,0,1,0)*q_saf(:,1) + &
                &VMEC_R(5,:,0,0,0)*0.0_dp,&
                &'check on DrtRq at par = 5')
            ! Drzeta
            Rq = 0.0_dp
            call arr_mult(VMEC_R,q_saf,Rq,[1,0,1])
            call write_out(1,n_r,Rq(5,:),'DrzRq at par = 5')
            call write_out(1,n_r,VMEC_R(5,:,1,0,1)*q_saf(:,0) + &
                &VMEC_R(5,:,1,0,0)*0.0_dp + &
                &VMEC_R(5,:,0,0,1)*q_saf(:,1) + &
                &VMEC_R(5,:,0,0,0)*0.0_dp,&
                &'check on DrzRq at par = 5')
            ! Dthetatheta
            Rq = 0.0_dp
            call arr_mult(VMEC_R,q_saf,Rq,[0,2,0])
            call write_out(1,n_r,Rq(5,:),'DttRq at par = 5')
            call write_out(1,n_r,VMEC_R(5,:,0,2,0)*q_saf(:,0) + &
                &2.0_dp * VMEC_R(5,:,0,1,0)*0.0_dp + &
                &VMEC_R(5,:,0,0,0)*0.0_dp,&
                &'check on DttRq at par = 5')
            ! Dtzeta
            Rq = 0.0_dp
            call arr_mult(VMEC_R,q_saf,Rq,[0,1,1])
            call write_out(1,n_r,Rq(5,:),'DtzRq at par = 5')
            call write_out(1,n_r,VMEC_R(5,:,0,1,1)*q_saf(:,0) + &
                &VMEC_R(5,:,0,1,0)*0.0_dp + &
                &VMEC_R(5,:,0,0,1)*0.0_dp + &
                &VMEC_R(5,:,0,0,0)*0.0_dp,&
                &'check on DtzRq at par = 5')
            ! Dzetazeta
            Rq = 0.0_dp
            call arr_mult(VMEC_R,q_saf,Rq,[0,0,2])
            call write_out(1,n_r,Rq(5,:),'DzzRq at par = 5')
            call write_out(1,n_r,VMEC_R(5,:,0,0,2)*q_saf(:,0) + &
                &2.0_dp * VMEC_R(5,:,0,0,1)*0.0_dp + &
                &VMEC_R(5,:,0,0,0)*0.0_dp,&
                &'check on DzzRq at par = 5')
            call lvl_ud(-1)
            
            ! higher order derive multiplied values
            call writo('higher order derive Rq')
            call lvl_ud(1)
            Rq = 0.0_dp
            call arr_mult(VMEC_R,q_saf,Rq,[1,2,3])
            call write_out(1,n_r,Rq(5,:),'DrttzzzRq at par = 5')
            call write_out(1,n_r,VMEC_R(5,:,1,2,3)*q_saf(:,0) + &
                &3.0_dp * VMEC_R(5,:,1,2,2)*0.0_dp + &
                &3.0_dp * VMEC_R(5,:,1,2,1)*0.0_dp + &
                &1.0_dp * VMEC_R(5,:,1,2,0)*0.0_dp + &
                &2.0_dp * VMEC_R(5,:,1,1,3)*0.0_dp + &
                &6.0_dp * VMEC_R(5,:,1,1,2)*0.0_dp + &
                &6.0_dp * VMEC_R(5,:,1,1,1)*0.0_dp + &
                &2.0_dp * VMEC_R(5,:,1,1,0)*0.0_dp + &
                &1.0_dp * VMEC_R(5,:,1,0,3)*0.0_dp + &
                &3.0_dp * VMEC_R(5,:,1,0,2)*0.0_dp + &
                &3.0_dp * VMEC_R(5,:,1,0,1)*0.0_dp + &
                &1.0_dp * VMEC_R(5,:,1,0,0)*0.0_dp + &
                &1.0_dp * VMEC_R(5,:,0,2,3)*q_saf(:,1) + &
                &3.0_dp * VMEC_R(5,:,0,2,2)*0.0_dp + &
                &3.0_dp * VMEC_R(5,:,0,2,1)*0.0_dp + &
                &1.0_dp * VMEC_R(5,:,0,2,0)*0.0_dp + &
                &2.0_dp * VMEC_R(5,:,0,1,3)*0.0_dp + &
                &6.0_dp * VMEC_R(5,:,0,1,2)*0.0_dp + &
                &6.0_dp * VMEC_R(5,:,0,1,1)*0.0_dp + &
                &2.0_dp * VMEC_R(5,:,0,1,0)*0.0_dp + &
                &1.0_dp * VMEC_R(5,:,0,0,3)*0.0_dp + &
                &3.0_dp * VMEC_R(5,:,0,0,2)*0.0_dp + &
                &3.0_dp * VMEC_R(5,:,0,0,1)*0.0_dp + &
                &1.0_dp * VMEC_R(5,:,0,0,0)*0.0_dp ,&
                &'check on DrttzzzRq at par = 5')
            call lvl_ud(-1)
        end subroutine
        
        subroutine case_1_1
            ! local variables
            real(dp) :: pq(n_r), pq_num(n_r)
            integer :: jd
            
            ! calculate equilibrium quantities
            call writo('calculate RZL')
            call lvl_ud(1)
            if (allocated(zeta)) deallocate(zeta)
            allocate(zeta(n_par,n_r)); zeta = 0.0_dp
            if (allocated(theta)) deallocate(theta)
            allocate(theta(n_par,n_r)); theta = 0.0_dp
            
            call init_eq
            
            zeta = 0.4*pi/2
            do jd = 1,n_r
                theta(:,jd) = eqd_mesh(n_par, 0.0_dp*pi, 3.0_dp*pi)
            end do
            call lvl_ud(-1)
            
            call writo('calculate q_saf and pres')
            call lvl_ud(1)
            call calc_flux_q
            call lvl_ud(-1)
            
            ! multiply
            call writo('multiply pres and q_saf')
            call lvl_ud(1)
            pq = 0.0_dp
            call arr_mult(pres,q_saf,pq,[0,0,0])
            pq_num = pres(:,0)*q_saf(:,0)
            call write_out(1,n_r,pq,'pq at par = 5')
            call write_out(1,n_r,pq_num,'check on pq at par = 5',comment=&
                &'max diff. = '//trim(r2str(maxval(abs(pq-pq_num)))))
            call lvl_ud(-1)
            
            ! derive multiplied values
            call writo('derive pq')
            call lvl_ud(1)
            ! Dr
            pq = 0.0_dp
            call arr_mult(pres,q_saf,pq,[1,0,0])
            pq_num = pres(:,1)*q_saf(:,0) + &
                &pres(:,0)*q_saf(:,1)
            call write_out(1,n_r,pq,'Drpq at par = 5')
            call write_out(1,n_r,pq_num,'check on Drpq at par = 5',comment=&
                &'max diff. = '//trim(r2str(maxval(abs(pq-pq_num)))))
            ! Dtheta
            pq = 0.0_dp
            call arr_mult(pres,q_saf,pq,[0,1,0])
            pq_num = 0.0_dp*q_saf(:,0) + &
                &pres(:,0)*0.0_dp
            call write_out(1,n_r,pq,'Dtpq at par = 5')
            call write_out(1,n_r,pq_num,'check on Dtpq at par = 5',comment=&
                &'max diff. = '//trim(r2str(maxval(abs(pq-pq_num)))))
            ! Dzeta
            pq = 0.0_dp
            call arr_mult(pres,q_saf,pq,[0,0,1])
            pq_num = 0.0_dp*q_saf(:,0) + &
                &pres(:,0)*0.0_dp
            call write_out(1,n_r,pq,'Dzpq at par = 5')
            call write_out(1,n_r,pq_num,'check on Dzpq at par = 5',comment=&
                &'max diff. = '//trim(r2str(maxval(abs(pq-pq_num)))))
            call lvl_ud(-1)
            
            ! double derive multiplied values
            call writo('double derive pq')
            call lvl_ud(1)
            ! Drr
            pq = 0.0_dp
            call arr_mult(pres,q_saf,pq,[2,0,0])
            call write_out(1,n_r,pq,'Drrpq at par = 5')
            call write_out(1,n_r,pres(:,2)*q_saf(:,0) + &
                &2.0_dp * pres(:,1)*q_saf(:,1) + &
                &pres(:,0)*q_saf(:,2),&
                &'check on Drrpq at par = 5')
            ! Drtheta
            pq = 0.0_dp
            call arr_mult(pres,q_saf,pq,[1,1,0])
            call write_out(1,n_r,pq,'Drrpq at par = 5')
            call write_out(1,n_r,0.0_dp*pq,&
                &'check on Drtpq at par = 5')
            ! Drzeta
            pq = 0.0_dp
            call arr_mult(pres,q_saf,pq,[1,0,1])
            call write_out(1,n_r,pq,'Drzpq at par = 5')
            call write_out(1,n_r,0.0_dp*pq,&
                &'check on Drzpq at par = 5')
            ! Dthetatheta
            pq = 0.0_dp
            call arr_mult(pres,q_saf,pq,[0,2,0])
            call write_out(1,n_r,pq,'Dttpq at par = 5')
            call write_out(1,n_r,0.0_dp*pq,&
                &'check on Dttpq at par = 5')
            ! Dtzeta
            pq = 0.0_dp
            call arr_mult(pres,q_saf,pq,[0,1,1])
            call write_out(1,n_r,pq,'Dtzpq at par = 5')
            call write_out(1,n_r,0.0_dp*pq,&
                &'check on Dtzpq at par = 5')
            ! Dzetazeta
            pq = 0.0_dp
            call arr_mult(pres,q_saf,pq,[0,0,2])
            call write_out(1,n_r,pq,'Dzzpq at par = 5')
            call write_out(1,n_r,0.0_dp*pq,&
                &'check on Dzzpq at par = 5')
            call lvl_ud(-1)
            
            ! higher order derive multiplied values
            call writo('higher order derive pq')
            call lvl_ud(1)
            pq = 0.0_dp
            call arr_mult(pres,q_saf,pq,[1,2,3])
            call write_out(1,n_r,pq,'Drttzzzpq at par = 5')
            call write_out(1,n_r,0.0_dp*pq ,&
                &'check on Drttzzzpq at par = 5')
            call lvl_ud(-1)
        end subroutine
    end subroutine
end module test
