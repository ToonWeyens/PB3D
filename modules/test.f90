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
        &test_ang_B, test_h2f_f2h, test_calc_norm_deriv, test_ext_var, test_B, &
        &test_fun_mult, test_norm_deriv, test_VMEC_jac

contains
    subroutine test_repack()
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

    subroutine test_write_out()
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

    subroutine test_mesh_cs()
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
    
    subroutine test_norm_deriv()                                                ! function version
        use utilities, only: norm_deriv
        use VMEC_vars, only: n_r
        
        ! local variables
        integer :: n_r_OLD
        integer :: n_max
        integer :: id
        
        call writo('test norm_deriv?')
        do
            if(yes_no(.false.)) then
                do
                    write(*,'(A)',advance='no') 'n_r ? ' 
                    read(*,*) n_max
                    if (n_max.lt.10) then
                        call writo('Choose a value larger than 10')
                        cycle
                    else
                        n_r_OLD = n_r
                        n_r = n_max
                        exit
                    end if
                end do
                
                call writo('FM values')
                call writo('Paused... press enter')
                read(*,*)
                call write_out(1,n_max,[(f(args(id),[0,0]),id=1,n_max)],&
                    &'f = cos x + 3*sin 2*y')
                call write_out(1,n_max,[(f100(args(id)),id=1,n_max)],&
                    &'Df = - sin x + 3*sin 2*y')
                call write_out(1,n_max,[(f100_num(args(id)),id=1,n_max)],&
                    &'Df num')
                call write_out(1,n_max,[(f100_diff(args(id)),id=1,n_max)],&
                    &'Df diff')
                
                call writo('HM values')
                call writo('Paused... press enter')
                read(*,*)
                call write_out(1,n_max-1,[(f(args_HM(id),[0,0]),id=2,n_max)],&
                    &'f = cos x + 3*sin 2*y')
                call write_out(1,n_max-1,[(f100(args_HM(id)),id=2,n_max)],&
                    &'Df = - sin x + 3*sin 2*y')
                call write_out(1,n_max-1,[(f100_num(args_HM(id)),id=2,n_max)],&
                    &'Df num')
                call write_out(1,n_max-1,[(f100_diff(args_HM(id)),id=2,n_max)],&
                    &'Df diff')
                
                n_r = n_r_OLD
            else
                exit
            end if
        end do
    contains
        function args(id)
            real(dp) :: args(3)
            integer, intent(in) :: id
            
            args = [1.0_dp*id,0.0_dp,0.0_dp]
        end function 
        function args_HM(id)
            real(dp) :: args_HM(3)
            integer, intent(in) :: id
            
            args_HM = [1.0_dp*id-0.5_dp,0.0_dp,0.0_dp]
        end function 
        function rad(pt)                                                        ! for pt = 1,n_max, returns 0,2*pi
            real(dp) :: rad
            real(dp), intent(in) :: pt(3)
            
            rad = 2*pi*(pt(1) - 1.0_dp)/(n_max-1.0_dp)
        end function rad
        function f(pt,deriv)
            real(dp) :: f
            real(dp), intent(in) :: pt(3)
            integer, intent(in) :: deriv(2)
            
            if (deriv(1).eq.0) then
                f = cos(rad(pt)) + 3*sin(2*rad(pt))
            else
                write(*,*) 'ERROR: only deriv(1) can be used'
            end if
        end function
        function f100(pt)
            real(dp) :: f100
            real(dp), intent(in) :: pt(3)
            
            f100 = - sin(rad(pt)) + 6*cos(2*rad(pt))
        end function
        function f100_num(pt)
            real(dp) :: f100_num
            real(dp), intent(in) :: pt(3)
            
            f100_num = norm_deriv(f,rad,pt,[0,0])
        end function
        function f100_diff(pt)
            real(dp) :: f100_diff
            real(dp), intent(in) :: pt(3)
            
            f100_diff = f100_num(pt) - f100(pt)
        end function
    end subroutine
    
    subroutine test_calc_norm_deriv()                                           ! array version
        use utilities, only: calc_norm_deriv
        
        ! local variables
        real(dp), allocatable :: varin(:), varout(:), varout_num(:)
        real(dp), allocatable :: x_FM(:), x_HM(:)
        integer :: id, n_r, tp
        real(dp), allocatable :: plotvar(:,:)
        
        do
            call writo('test calc_norm_deriv?')
            if(yes_no(.false.)) then
                write(*,'(A)',advance='no') 'n_r ? ' 
                read(*,*) n_r
                write(*,*) 'type ?'
                write(*,'(A)',advance='no') '(1: (FM_i,FM_o), 2: (FM_i,HM_o), &
                    &3: (HM_i,FM_o), 4: (HM_i,HM_o)  = ' 
                read(*,*) tp
                
                if (allocated(x_FM)) deallocate(x_FM)
                if (allocated(x_HM)) deallocate(x_HM)
                if (allocated(varin)) deallocate(varin)
                if (allocated(varout)) deallocate(varout)
                if (allocated(varout_num)) deallocate(varout_num)
                if (allocated(plotvar)) deallocate(plotvar)
                allocate(x_FM(n_r))
                allocate(x_HM(n_r))
                allocate(varin(n_r))
                allocate(varout(n_r))
                allocate(varout_num(n_r))
                allocate(plotvar(2, n_r))
                
                ! x = [0...1]
                x_FM = [(float(id)/(n_r-1), id=0,n_r-1)]
                x_HM = x_FM - 0.5/(n_r-1)
                
                ! input: sin(2*pi * x) + 0.5*cos(2*pi * 0.3*x)
                if (tp.eq.1 .or. tp.eq.2) then                                  ! FM output
                    varin = sin(2*pi*x_FM) + 0.5* cos(2*pi*0.3*x_FM)
                else                                                            ! HM output
                    varin = sin(2*pi*x_HM) + 0.5* cos(2*pi*0.3*x_HM)
                    varin(1) = 0.0_dp
                end if
                
                ! output = 2*pi * cos(2*pi * x)
                if (tp.eq.1 .or. tp.eq.3) then                                  ! FM output
                    varout = 2*pi * (cos(2*pi*x_FM) - 0.15*sin(2*pi*0.3*x_FM))
                else                                                            ! HM output
                    varout = 2*pi * (cos(2*pi*x_HM) - 0.15*sin(2*pi*0.3*x_HM))
                    varout(1) = 0.0_dp
                end if
                
                select case (tp)
                    case (1)                                                    ! FM input, FM output
                        varout_num = calc_norm_deriv(varin,&
                            &1.0_dp*float(n_r-1),.true.,.true.)
                    case (2)                                                    ! FM input, HM output
                        varout_num = calc_norm_deriv(varin,&
                            &1.0_dp*float(n_r-1),.true.,.false.)
                        varout(1) = 0.0_dp
                    case (3)                                                    ! HM input, FM output
                        varout_num = calc_norm_deriv(varin,&
                            &1.0_dp*float(n_r-1),.false.,.true.)
                    case (4)                                                    ! HM input, HM output
                        varout_num = calc_norm_deriv(varin,&
                            &1.0_dp*float(n_r-1),.false.,.false.)
                        varout(1) = 0.0_dp
                end select
                
                if (tp.eq.1 .or. tp.eq.2) then                                  ! FM input
                    plotvar(1,:) = x_FM
                    plotvar(2,:) = varin
                else
                    plotvar(1,:) = x_HM
                    plotvar(2,:) = varin
                end if
                call write_out(2,n_r,plotvar,'varin')
                
                if (tp.eq.1 .or. tp.eq.3) then                                  ! FM output
                    plotvar(1,:) = x_FM
                    plotvar(2,:) = varout
                else
                    plotvar(1,:) = x_HM
                    plotvar(2,:) = varout
                end if
                call write_out(2,n_r,plotvar,'varout')
                
                if (tp.eq.1 .or. tp.eq.3) then                                  ! FM output
                    plotvar(1,:) = x_FM
                    plotvar(2,:) = varout_num
                else
                    plotvar(1,:) = x_HM
                    plotvar(2,:) = varout_num
                end if
                call write_out(2,n_r,plotvar,'varout numerical')
                
                if (tp.eq.1 .or. tp.eq.3) then                                  ! FM output
                    plotvar(1,:) = x_FM
                    plotvar(2,:) = abs(varout_num-varout)
                else
                    plotvar(1,:) = x_HM
                    plotvar(2,:) = abs(varout_num-varout)
                end if
                call write_out(2,n_r,plotvar,'absolute difference')
            else
                exit
            end if
        end do
    end subroutine
    
    subroutine test_ext_var()
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
    
    subroutine test_h2f_f2h()
        use utilities, only: h2f, f2h
        
        ! local variables
        integer :: id, kd
        integer :: n_max, step, n_r
        integer, parameter :: length = 50
        real(dp), allocatable :: varin(:), varout(:), varoutout(:), vardiff(:)
        real(dp) :: maxerr(length), averr(length), num_points(length)
        real(dp) :: plotvar(2,length)
        logical :: ind_plot, log_plot
        
        call writo('test h2f, f2h?')
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
                
                varout = f2h(varin)
                varoutout = h2f(varout)
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
    
    subroutine test_VMEC_jac()
        use eq_ops, only: calc_eq
        use eq_vars, only: n_par, max_par, min_par, theta_H, zeta_H
        use VMEC_vars, only: n_r
        use metric_ops, only: jac_V_H, VMEC_jac
        
        ! local variables
        real(dp) :: jac_V_H_alt2(n_par, n_r)                                    ! Jacobian from functions (HM)
        real(dp) :: jac_V_H_der(n_par, n_r)                                     ! angular deriv. of Jacobian from functions (HM)
        real(dp) :: jac_V_H_diff(n_r)                                           ! difference over all parallel points
        integer :: jd,kd
        
        call writo('test VMEC_jac?')
        if(yes_no(.false.)) then
            ! calculate equilibrium with some value for alpha
            call writo('Calculate equilibrium')
            call lvl_ud(1)
            call calc_eq(pi/2)
            call lvl_ud(-1)
            
            ! calculate jacobian using VMEC_jac
            call writo('Calculating jacobians')
            jac_V_H_alt2 = 0.0_dp
            par: do jd = 1, n_par                                               ! parallel: along the magnetic field line
                perp: do kd = 2, n_r                                            ! perpendicular: normal to the flux surfaces
                    jac_V_H_alt2(jd,kd) = VMEC_jac([kd-0.5_dp,theta_H(jd,kd),&
                        &zeta_H(jd,kd)],[0,0])
                end do perp
            end do par
            
            ! plot
            call writo('Plotting jacobians')
            call lvl_ud(1)
            call writo('individual plots?')
            if(yes_no(.false.)) then
                do jd = 1,n_par
                    call write_out(1,n_r,jac_V_H(jd,:),'jac_V_H at par = '//&
                        &trim(i2str(jd)))
                    call write_out(1,n_r,jac_V_H_alt2(jd,:),'jac_V_H_alt2 at &
                        &par = '//trim(i2str(jd)))
                    call write_out(1,n_r,abs(jac_V_H_alt2(jd,:)-jac_V_H(jd,:)),&
                        &'jac_V_H_alt2 at par = '//trim(i2str(jd)))
                end do
            end if
            do kd = 1,n_r
                jac_V_H_diff(kd) = maxval(abs(jac_V_H_alt2(:,kd)-jac_V_H(:,kd)))
            end do
            call write_out(1,n_r,jac_V_H_diff,'maximum difference between &
                &jac_V_H')
            call lvl_ud(-1)
            
            ! angular derivatives
            call writo('Calculating angular derivatives')
            call lvl_ud(1)
            perp2: do kd = 2, n_r
                call writo('Poloidal dependency at r = '//trim(i2str(kd)))
                parpol: do jd = 1, n_par
                    jac_V_H_alt2(jd,kd) = VMEC_jac([kd-0.5_dp,(max_par-min_par)&
                        &*jd/n_par,0.4*pi],[0,0])
                    jac_V_H_der(jd,kd) = VMEC_jac([kd-0.5_dp,(max_par-min_par)&
                        &*jd/n_par,0.4*pi],[1,0])
                end do parpol
                call write_out(1,n_par,jac_V_H_alt2(:,kd),'jacobian at r = '&
                    &//trim(i2str(kd)))
                call write_out(1,n_par,jac_V_H_der(:,kd),'pol. deriv at r = '&
                    &//trim(i2str(kd)))
                
                call writo('Toroidal dependency at r = '//trim(i2str(kd)))
                partor: do jd = 1, n_par
                    jac_V_H_alt2(jd,kd) = VMEC_jac([kd-0.5_dp,0.4*pi,&
                        &(max_par-min_par)*jd/n_par],[0,0])
                    jac_V_H_der(jd,kd) = VMEC_jac([kd-0.5_dp,0.4*pi,&
                        &(max_par-min_par)*jd/n_par],[0,1])
                end do partor
                call write_out(1,n_par,jac_V_H_alt2(:,kd),'jacobian at r = '&
                    &//trim(i2str(kd)))
                call write_out(1,n_par,jac_V_H_der(:,kd),'tor. deriv at r = '&
                    &//trim(i2str(kd)))
            end do perp2
            call lvl_ud(-1)
        end if
    end subroutine

    subroutine test_metric_transf()
        use metric_ops, only: metric_C, metric_C2V, metric_V, metric_V2F, &
            &metric_F, &
            &C2V_up, C2V_dn, C2V_up_H, C2V_dn_H, V2F_up, V2F_dn, V2F_up_H, &
            &V2F_dn_H, jac_V, g_V, h_V, jac_F, g_F, h_F, jac_V_H, g_V_H, &
            &h_V_H, jac_F_H, g_F_H, h_F_H
        use eq_vars, only: eqd_mesh, calc_RZL, calc_flux_q, &
            &n_par, min_par, max_par, theta, zeta, theta_H, zeta_H
        use VMEC_vars, only: n_r
        use utilities, only: mat_mult, det
        
        ! local variables
        integer :: kd
        
        call writo('test metric_transf?')
        if(yes_no(.false.)) then
            ! initalize theta and zeta. No need for field-line following theta
            if (allocated(zeta)) deallocate(zeta)
            if (allocated(theta)) deallocate(theta)
            if (allocated(zeta_H)) deallocate(zeta_H)
            if (allocated(theta_H)) deallocate(theta_H)
            allocate(zeta(n_par, n_r)); zeta = 0.0_dp
            allocate(theta(n_par, n_r)); theta = 0.0_dp
            allocate(zeta_H(n_par, n_r)); zeta_H = 0.0_dp
            allocate(theta_H(n_par, n_r)); theta_H = 0.0_dp
            do kd = 1,n_r
                zeta(:,kd) = eqd_mesh(n_par, min_par, max_par)
                zeta_H(:,kd) = eqd_mesh(n_par, min_par, max_par)
            end do
            theta = zeta_H
            theta_H = zeta_H
            
            ! calculate the cylindrical variables R, Z and lambda and derivatives
            call calc_RZL
            
            ! calculate flux quantities
            call calc_flux_q
            
            ! calculate the metrics in the cylindrical coordinate system
            call metric_C
            
            ! calculate the transformation matrix C(ylindrical) -> V(mec)
            call metric_C2V
            
            ! calculate  the  metric  factors in the VMEC coordinate system 
            call metric_V
            
            ! calculate the transformation matrix V(mec) -> F(lux)
            call metric_V2F
            
            ! calculate  the  metric  factors in the flux coordinate system 
            call metric_F
            
            call writo('test whether we have C2V_dn*C2V_up^T = 1')
            call lvl_ud(1)
            
            call writo('Cylindrical -> VMEC:')
            call lvl_ud(1)
            call writo('Half Mesh:')
            call lvl_ud(1)
            call writo('Paused... press enter')
            read(*,*)
            call lvl_ud(-1)
            call transf_inv_transf(C2V_up_H, C2V_dn_H)
            call writo('Full Mesh:')
            call lvl_ud(1)
            call writo('Paused... press enter')
            read(*,*)
            call lvl_ud(-1)
            call transf_inv_transf(C2V_up, C2V_dn)
            call lvl_ud(-1)
            call writo('VMEC -> flux:')
            call lvl_ud(1)
            call writo('Half Mesh:')
            call lvl_ud(1)
            call writo('Paused... press enter')
            read(*,*)
            call transf_inv_transf(V2F_up_H, V2F_dn_H)
            call lvl_ud(-1)
            call writo('Full Mesh:')
            call lvl_ud(1)
            call writo('Paused... press enter')
            read(*,*)
            call transf_inv_transf(V2F_up, V2F_dn)
            call lvl_ud(-1)
            call lvl_ud(-1)
            
            call lvl_ud(-1)
            
            call writo('Test whether the determinants of g_C, h_C are agreee &
                &with jac_V')
            call lvl_ud(1)
            
            call writo('VMEC:')
            call lvl_ud(1)
            call writo('Half Mesh:')
            call lvl_ud(1)
            call writo('Paused... press enter')
            read(*,*)
            call check_dets(jac_V_H,g_V_H,h_V_H)
            call lvl_ud(-1)
            call writo('Full Mesh:')
            call lvl_ud(1)
            call writo('Paused... press enter')
            read(*,*)
            call check_dets(jac_V,g_V,h_V)
            call lvl_ud(-1)
            call lvl_ud(-1)
            call writo('flux:')
            call lvl_ud(1)
            call writo('Half Mesh:')
            call lvl_ud(1)
            call writo('Paused... press enter')
            read(*,*)
            call check_dets(jac_F_H,g_F_H,h_F_H)
            call lvl_ud(-1)
            call writo('Full Mesh:')
            call lvl_ud(1)
            call writo('Paused... press enter')
            read(*,*)
            call check_dets(jac_F,g_F,h_F)
            call lvl_ud(-1)
            call lvl_ud(-1)
            
            call lvl_ud(-1)
            
            call writo('Paused... press enter')
            read(*,*)
        end if
    contains
        ! calculate T_dn * T_up^T and check if it's equal to 1
        subroutine transf_inv_transf(T_up,T_dn)
            ! input / output
            real(dp) :: T_up(3,3,n_par,n_r), T_dn(3,3,n_par,n_r) 
            
            ! local variables
            real(dp) :: T_mult(3,3)
            real(dp) :: u_mat(3,3), diff_mat(3,3)                               ! unity 3x3 matrix, difference matrix with unity_mat
            real(dp) :: diff_max(n_r)                                           ! maximum deviation of unity matrix for all radial points
            integer :: id
            
            u_mat = 0.0_dp
            u_mat(1,1) = 1.0_dp; u_mat(2,2) = 1.0_dp; u_mat(3,3) = 1.0_dp
            diff_max = 0.0_dp
            do kd = 2,n_r                                                       ! the magnetic axis gives infinities
                do id = 1,n_par
                    T_mult = mat_mult(T_up(:,:,id,kd),&
                        &transpose(T_dn(:,:,id,kd)))
                    diff_mat = T_mult - u_mat
                    if (maxval(abs(diff_mat)).gt.diff_max(kd)) then
                        diff_max(kd) = log10(max(maxval(abs(diff_mat)),&
                            &1.0E-20_dp))
                    end if
                end do
            end do
            call write_out(1,n_r,diff_max,'maximum difference from &
                &unit matrix [log]', comment=' for different radial positions')
        end subroutine
        
        ! check whether jacobians calculated in the routines in metric_ops and
        ! the correct determinants are in agreement
        subroutine check_dets(jac,g,h)
            ! input / output
            real(dp) :: jac(n_par,n_r)
            real(dp) :: g(3,3,n_par,n_r), h(3,3,n_par,n_r)
            
            ! local variables
            real(dp) :: g_J, h_J                                                ! jacobian calculated from g and h
            real(dp) :: diff_J_max(n_r,2)                                       ! maximum of difference between g_J, h_J and jac_V
            integer :: id, kd
            !real(dp) :: j_calc(n_par,n_r,2)
            
            diff_J_max = 0.0_dp
            do kd = 1,n_r
                do id = 1,n_par
                    g_J = sqrt(abs(det(3,g(:,:,id,kd))))
                    !j_calc(id,kd,1) = g_J
                    h_J = 1.0_dp/sqrt(abs(det(3,h(:,:,id,kd))))
                    !j_calc(id,kd,2) = h_J
                    if (abs(abs(g_J)-abs(jac(id,kd))).gt.diff_J_max(kd,1)) then
                        diff_J_max(kd,1) = log10(max(abs(abs(g_J)-abs(jac(id,kd))),&
                            &1.0E-20_dp))
                    end if
                    if (abs(abs(h_J)-abs(jac(id,kd))).gt.diff_J_max(kd,2)) then
                        diff_J_max(kd,2) = log10(max(abs(abs(h_J)-abs(jac(id,kd))),&
                            &1.0E-20_dp))
                    end if
                end do
            end do
            
            call write_out(1,n_r,diff_J_max(:,1),'maximum difference &
                &between sqrt(det(g_ij)) and jac [log]', &
                &comment=' for different radial positions')
            call write_out(1,n_r,diff_J_max(:,2),'maximum difference &
                &between sqrt(det(h^ij)) and jac [log]', &
                &comment=' for different radial positions')
        end subroutine
    end subroutine

    subroutine test_B()
        use eq_ops, only: calc_eq
        use eq_vars, only: n_par, flux_p, flux_p_H, q_saf, q_saf_H, lam_H
        use metric_ops, only: jac_V, jac_V_H, g_V, g_V_H, jac_F, jac_F_H
        use VMEC_vars, only: n_r
        use utilities, only: h2f
        use B_vars, only: B_V_sub, B_V_sub_H, B_V, B_V_H
        
        ! local variables
        real(dp) :: B_V_sub_alt(n_par,n_r,3)
        real(dp) :: B_V_sub_alt2(n_par,n_r,3)
        real(dp) :: B_V_sub_H_alt(n_par,n_r,3)
        real(dp) :: B_V_sub_H_alt2(n_par,n_r,3)
        real(dp) :: B_V_alt(n_par,n_r)
        real(dp) :: B_V_H_alt(n_par,n_r)
        real(dp) :: diff(n_r)
        integer :: id, jd, kd
        real(dp) :: lam(n_par,n_r,4)
        
        call writo('test magnetic fields?')
        if(yes_no(.false.)) then
            call lvl_ud(1)
            
            ! calculate equilibrium with some value for alpha
            call writo('Calculate equilibrium')
            call lvl_ud(1)
            call calc_eq(pi/2)
            call lvl_ud(-1)
            
            call writo('Testing whether the magnetic fields from VMEC are &
                &in agreement with B_i = nabla alpha x nabla psi')
            call lvl_ud(1)
            
            ! calculate lam from lam_H
            do id = 1,n_par
                do jd = 1,4
                    lam(id,:,jd) = h2f(lam_H(id,:,jd))
                end do
            end do
            
            ! calculate B_i alternatively
            B_V_sub_alt = 0.0_dp
            B_V_sub_H_alt = 0.0_dp
            perp: do kd = 1,n_r                                                 ! avoid problems at magn. axis.
                par: do id = 1,n_par
                    comp: do jd = 1,3
                        B_V_sub_alt(id,kd,jd) = &
                            &flux_p(kd,2)/(2*pi*jac_V(id,kd)) * &
                            &(q_saf(kd,1)*(1+lam(id,kd,3)) * g_V(3,jd,id,kd) + &
                            &(1-q_saf(kd,1)*lam(id,kd,4)) * g_V(2,jd,id,kd))
                        B_V_sub_alt2(id,kd,jd) = &
                            &1/jac_F(id,kd) * &
                            &(q_saf(kd,1) * g_V(3,jd,id,kd) + &
                            &(1-q_saf(kd,1)*lam(id,kd,4))/(1+lam(id,kd,3))* &
                            &g_V(2,jd,id,kd))
                        B_V_sub_H_alt(id,kd,jd) = &
                            &flux_p_H(kd,2)/(2*pi*jac_V_H(id,kd)) * &
                            &(q_saf_H(kd,1)*(1+lam_H(id,kd,3)) &
                            &* g_V_H(3,jd,id,kd) + (1-q_saf_H(kd,1)*&
                            &lam_H(id,kd,4)) * g_V_H(2,jd,id,kd))
                        B_V_sub_H_alt2(id,kd,jd) = &
                            &1/jac_F_H(id,kd) * &
                            &(q_saf_H(kd,1) * g_V_H(3,jd,id,kd) + &
                            &(1-q_saf_H(kd,1)*lam_H(id,kd,4))/&
                            &(1+lam_H(id,kd,3))* g_V_H(2,jd,id,kd))
                    end do comp
                end do par
            end do perp
            
            ! calculate B alternatively
            B_V_alt = 0.0_dp
            B_V_H_alt = 0.0_dp
            perp2: do kd = 1,n_r
                par2: do id = 1,n_par
                    B_V_alt(id,kd) = sqrt(1/(jac_F(id,kd)) * &
                        &((1-q_saf(kd,1)*lam(id,kd,4))/(1+lam(id,kd,3)) &
                        &* B_V_sub(id,kd,1,2) + q_saf(kd,1)*B_V_sub(id,kd,1,3)))
                    B_V_H_alt(id,kd) = sqrt(1/(jac_F_H(id,kd)) * &
                        &((1-q_saf_H(kd,1)*lam_H(id,kd,4))/(1+lam_H(id,kd,3)) &
                        &* B_V_sub_H(id,kd,1,2) &
                        &+ q_saf_H(kd,1)*B_V_sub_H(id,kd,1,3)))
                end do par2
            end do perp2
            
            call writo('FULL MESH')
            call lvl_ud(1)
            diff = 0.0_dp
            call writo('Paused... press enter')
            read(*,*)
            
            ! FM magnitude
            do kd = 1, n_r
                diff(kd) = log10(min(max(maxval(abs(B_V(:,kd)-B_V_alt(:,kd)))&
                    &,1.0E-20_dp),1E10_dp))
            end do
            call write_out(1,n_r,diff,'test whether sqrt(B_theta/J)_F = B &
                &[log]')
            
            ! FM r component
            do kd = 1, n_r
                diff(kd) = log10(min(max(maxval(abs(B_V_sub(:,kd,1,1)-&
                    &B_V_sub_alt(:,kd,1))),1.0E-20_dp),1E10_dp))
            end do
            call write_out(1,n_r,diff,'max diff. between standard and &
                &alternative treatment, r component [log]')
            do kd = 1, n_r
                diff(kd) = log10(min(max(maxval(abs(B_V_sub(:,kd,1,1)-&
                    &B_V_sub_alt2(:,kd,1))),1.0E-20_dp),1E10_dp))
            end do
            call write_out(1,n_r,diff,'max diff. between standard and &
                &alternative treatment 2, r component [log]')
            
            ! FM theta component
            do kd = 1, n_r
                diff(kd) = log10(min(max(maxval(abs(B_V_sub(:,kd,1,2)-&
                    &B_V_sub_alt(:,kd,2))),1.0E-20_dp),1E10_dp))
            end do
            call write_out(1,n_r,diff,'max diff. between standard and &
                &alternative treatment, theta component [log]')
            do kd = 1, n_r
                diff(kd) = log10(min(max(maxval(abs(B_V_sub(:,kd,1,2)-&
                    &B_V_sub_alt2(:,kd,2))),1.0E-20_dp),1E10_dp))
            end do
            call write_out(1,n_r,diff,'max diff. between standard and &
                &alternative treatment 2, theta component [log]')
                
            ! FM zeta component
            do kd = 1, n_r
                diff(kd) = log10(min(max(maxval(abs(B_V_sub(:,kd,1,3)-&
                    &B_V_sub_alt(:,kd,3))),1.0E-20_dp),1E10_dp))
            end do
            call write_out(1,n_r,diff,'max diff. between standard and &
                &alternative treatment, zeta component [log]')
            do kd = 1, n_r
                diff(kd) = log10(min(max(maxval(abs(B_V_sub(:,kd,1,3)-&
                    &B_V_sub_alt2(:,kd,3))),1.0E-20_dp),1E10_dp))
            end do
            call write_out(1,n_r,diff,'max diff. between standard and &
                &alternative treatment 2, zeta component [log]')
            
            call lvl_ud(-1)
            
            call writo('HALF MESH')
            call lvl_ud(1)
            call writo('Paused... press enter')
            read(*,*)
            
            ! HM magnitude
            do kd = 1, n_r
                diff(kd) = log10(min(max(maxval(abs(B_V_H(:,kd)-&
                    &B_V_H_alt(:,kd))),1.0E-20_dp),1E10_dp))
            end do
            call write_out(1,n_r-1,diff(2:n_r),'test whether &
                &sqrt(B_theta/J)_F = B [log]')
            
            ! HM r component
            do kd = 1, n_r
                diff(kd) = log10(max(maxval(abs(B_V_sub_H(:,kd,1,1)-&
                    &B_V_sub_H_alt(:,kd,1))),1.0E-20_dp))
            end do
            call write_out(1,n_r-1,diff(2:n_r),'max diff. between standard and &
                &alternative treatment, r component [log]')
            do kd = 1, n_r
                diff(kd) = log10(max(maxval(abs(B_V_sub_H(:,kd,1,1)-&
                    &B_V_sub_H_alt2(:,kd,1))),1.0E-20_dp))
            end do
            call write_out(1,n_r-1,diff(2:n_r),'max diff. between standard and &
                &alternative treatment 2, r component [log]')
            
            ! HM theta component
            do kd = 1, n_r
                diff(kd) = log10(max(maxval(abs(B_V_sub_H(:,kd,1,2)-&
                    &B_V_sub_H_alt(:,kd,2))),1.0E-20_dp))
            end do
            call write_out(1,n_r-1,diff(2:n_r),'max diff. between standard and &
                &alternative treatment, theta component [log]')
            do kd = 1, n_r
                diff(kd) = log10(max(maxval(abs(B_V_sub_H(:,kd,1,2)-&
                    &B_V_sub_H_alt2(:,kd,2))),1.0E-20_dp))
            end do
            call write_out(1,n_r-1,diff(2:n_r),'max diff. between standard and &
                &alternative treatment 2, theta component [log]')
            
            ! HM zeta component
            do kd = 1, n_r
                diff(kd) = log10(max(maxval(abs(B_V_sub_H(:,kd,1,3)-&
                    &B_V_sub_H_alt(:,kd,3))),1.0E-20_dp))
            end do
            call write_out(1,n_r-1,diff(2:n_r),'max diff. between standard and &
                &alternative treatment, zeta component [log]')
            do kd = 1, n_r
                diff(kd) = log10(max(maxval(abs(B_V_sub_H(:,kd,1,3)-&
                    &B_V_sub_H_alt2(:,kd,3))),1.0E-20_dp))
            end do
            call write_out(1,n_r-1,diff(2:n_r),'max diff. between standard and &
                &alternative treatment 2, zeta component [log]')
            
            call lvl_ud(-1)
            
            call lvl_ud(-1)
            call lvl_ud(-1)
        end if
    end subroutine

    subroutine test_ang_B()
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
    
    subroutine test_fun_mult()
        use utilities, only: fun_mult
        use fourier_ops, only: f2r
        use eq_vars, only: VMEC_R, vmec_z
        use VMEC_vars, only: ntor
        
        ! local variables also used in subfunctinos
        real(dp) :: rad, theta, zeta
        integer :: n_pt                                                         ! max # points to plot
        integer :: id                                                           ! counter
        
        call writo('test fun_mult?')
        if(yes_no(.false.)) then
            call writo('Testing whether fun_mult correctly multiplies two &
                &functions, taking into account derivatives')
            call lvl_ud(1)
            
            ! read user-input
            write(*,'(A)',advance='no') 'r = ? ' 
            read(*,*) rad
            write(*,'(A)',advance='no') 'theta = ? ' 
            read(*,*) theta
            if (ntor.gt.0) then
                write(*,'(A)',advance='no') 'zeta = ? ' 
                read(*,*) zeta
            end if
            
            ! plot R
            n_pt = 50
            
            ! R
            call write_out(1,n_pt,[(R([rad,id*2*pi/n_pt,zeta]),&
                &id=1,n_pt)],'R at r = '//trim(r2strt(rad))//' and zeta = '//&
                &trim(r2strt(zeta)))
            ! Z
            call write_out(1,n_pt,[(Z([rad,id*2*pi/n_pt,zeta]),&
                &id=1,n_pt)],'Z at r = '//trim(r2strt(rad))//' and zeta = '//&
                &trim(r2strt(zeta)))
            ! R Z
            call write_out(1,n_pt,[(RZ([rad,id*2*pi/n_pt,zeta]),&
                &id=1,n_pt)],'R Z at r = '//trim(r2strt(rad))//&
                &' and zeta = '//trim(r2strt(zeta)))
            ! R Z_alt
            call write_out(1,n_pt,[(RZ_alt([rad,id*2*pi/n_pt,zeta]),&
                &id=1,n_pt)],'R Z alt at r = '//trim(r2strt(rad))//&
                &' and zeta = '//trim(r2strt(zeta)))
            ! R Z_diff
            call write_out(1,n_pt,[(RZ_diff([rad,id*2*pi/n_pt,zeta]),&
                &id=1,n_pt)],'R Z alt diff r = '//trim(r2strt(rad))//&
                &' and zeta = '//trim(r2strt(zeta)))
            ! R Z 30
            call write_out(1,n_pt,[(RZ30([rad,id*2*pi/n_pt,zeta]),&
                &id=1,n_pt)],'R Z 30 at r = '//trim(r2strt(rad))//&
                &' and zeta = '//trim(r2strt(zeta)))
            ! R Z_alt 30
            call write_out(1,n_pt,[(RZ30_alt([rad,id*2*pi/n_pt,zeta]),&
                &id=1,n_pt)],'R Z 30 alt at r = '//trim(r2strt(rad))//&
                &' and zeta = '//trim(r2strt(zeta)))
            ! R Z_diff 30
            call write_out(1,n_pt,[(RZ30_diff([rad,id*2*pi/n_pt,zeta]),&
                &id=1,n_pt)],'R Z 30 diff at r = '//trim(r2strt(rad))//&
                &' and zeta = '//trim(r2strt(zeta)))
                
            if (ntor.gt.0) then
                ! R Z 11
                call write_out(1,n_pt,[(RZ11([rad,id*2*pi/n_pt,zeta]),&
                    &id=1,n_pt)],'R Z 11 at r = '//trim(r2strt(rad))//&
                    &' and zeta = '//trim(r2strt(zeta)))
                ! R Z 11 alt
                call write_out(1,n_pt,[(RZ11_alt([rad,id*2*pi/n_pt,zeta]),&
                    &id=1,n_pt)],'R Z 11 alt at r = '//trim(r2strt(rad))//&
                    &' and zeta = '//trim(r2strt(zeta)))
                ! R Z 11 diff
                call write_out(1,n_pt,[(RZ11_diff([rad,id*2*pi/n_pt,zeta]),&
                    &id=1,n_pt)],'R Z 11 diff at r = '//trim(r2strt(rad))//&
                    &' and zeta = '//trim(r2strt(zeta)))
            end if
            
            call lvl_ud(-1)
        end if
    contains
        function R(pt)
            real(dp) :: R
            real(dp) :: pt(3)
            
            R = VMEC_R(pt,[0,0])
        end function
        function Z(pt)
            real(dp) :: Z
            real(dp) :: pt(3)
            
            Z = VMEC_Z(pt,[0,0])
        end function
        function RZ(pt)
            real(dp) :: RZ
            real(dp) :: pt(3)
            
            RZ = fun_mult(VMEC_R,VMEC_Z,pt,[0,0])
        end function
        function RZ_alt(pt)
            real(dp) :: RZ_alt
            real(dp) :: pt(3)
            
            RZ_alt = VMEC_R(pt,[0,0])*VMEC_Z(pt,[0,0])
        end function
        function RZ_diff(pt)
            real(dp) :: RZ_diff
            real(dp) :: pt(3)
            
            RZ_diff = VMEC_R(pt,[0,0])*VMEC_Z(pt,[0,0]) - RZ(pt)
        end function
        function RZ30(pt)
            real(dp) :: RZ30
            real(dp) :: pt(3)
            
            RZ30 = fun_mult(VMEC_R,VMEC_Z,pt,[3,0])
        end function
        function RZ30_alt(pt)
            real(dp) :: RZ30_alt
            real(dp) :: pt(3)
            
            RZ30_alt = VMEC_R(pt,[3,0])*VMEC_Z(pt,[0,0]) + 3*VMEC_R(pt,[2,0])&
                &*VMEC_Z(pt,[1,0]) + 3*VMEC_R(pt,[1,0])*VMEC_Z(pt,[2,0]) + &
                &VMEC_R(pt,[0,0])*VMEC_Z(pt,[3,0])
        end function
        function RZ30_diff(pt)
            real(dp) :: RZ30_diff
            real(dp) :: pt(3)
            
            RZ30_diff = RZ30_alt(pt) - RZ30(pt)
        end function
        function R01(pt)
            real(dp) :: R01
            real(dp) :: pt(3)
            
            R01 = VMEC_R(pt,[0,1])
        end function
        function R10(pt)
            real(dp) :: R10
            real(dp) :: pt(3)
            
            R10 = VMEC_R(pt,[1,0])
        end function
        function R11(pt)
            real(dp) :: R11
            real(dp) :: pt(3)
            
            R11 = VMEC_R(pt,[1,1])
        end function
        function Z01(pt)
            real(dp) :: Z01
            real(dp) :: pt(3)
            
            Z01 = VMEC_Z(pt,[0,1])
        end function
        function Z10(pt)
            real(dp) :: Z10
            real(dp) :: pt(3)
            
            Z10 = VMEC_Z(pt,[1,0])
        end function
        function Z11(pt)
            real(dp) :: Z11
            real(dp) :: pt(3)
            
            Z11 = VMEC_Z(pt,[1,1])
        end function
        function RZ11(pt)
            real(dp) :: RZ11
            real(dp) :: pt(3)
            
            RZ11 = fun_mult(VMEC_R,VMEC_Z,pt,[1,1])
        end function
        function RZ11_alt(pt)
            real(dp) :: RZ11_alt
            real(dp) :: pt(3)
            
            RZ11_alt = R11(pt)*Z(pt) + R01(pt)*Z10(pt) + R10(pt)*Z01(pt) &
                &+ R(pt)*Z11(pt)
        end function
        function RZ11_diff(pt)
            real(dp) :: RZ11_diff
            real(dp) :: pt(3)
            
            RZ11_diff = RZ11_alt(pt) - RZ11(pt)
        end function
    end subroutine
end module test
