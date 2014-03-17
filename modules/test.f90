! Test some routines and functions
module test
    use num_vars, only: dp, max_str_ln, pi
    use output_ops, only: writo, lvl_ud, print_ar_1, print_ar_2, write_out, &
        &lvl
    use time, only: start_time, stop_time
    use str_ops, only: strh2l, i2str, r2str, r2strt
    use var_ops, only: mat_mult, det

    implicit none
    private
    public test_repack, test_write_out, test_mesh_cs, test_metric_C2V, &
        &test_theta_B

contains
    subroutine test_repack()
        ! VMEC variable has structure (1:mnmax, 1:n_r)
        ! output variable should have (1:n_r, 0:mpol-1, -ntor:ntor)
        use fourier_ops, only: repack
        integer :: n_rB, mpolB, ntorB, mnmaxB
        real(dp), allocatable :: xmB(:), xnB(:)
        real(dp), allocatable :: varinB(:,:)
        real(dp), allocatable :: varoutB(:,:,:)
            
        if(test_this('repack')) then
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
            
        if(test_this('write_out')) then
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
            if(test_this('mesh_cs')) then
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

    subroutine test_metric_C2V()
        use metric_ops, only: metric_C, metric_C2V, metric_V, &
            &C2V_up, C2V_dn, jac_V, g_V, h_V
        use eq_vars, only: eqd_mesh, calc_RZl, n_zeta, &
            &min_zeta, max_zeta, theta, zeta
        use VMEC_vars, only: n_r
        
        real(dp) :: C2V_mult(3,3)
        real(dp) :: u_mat(3,3), diff_mat(3,3)                                   ! unity 3x3 matrix, difference matrix with unity_mat
        integer :: max_index(3)                                                 ! index of maximum difference
        real(dp) :: diff_max                                                    ! maximum deviation of unity matrix
        real(dp) :: g_J, h_J                                                    ! jacobian calculated from g and h
        integer :: max_J_index(2,3)                                             ! index of maximum difference 
        real(dp) :: diff_J_max(2)                                               ! maximum of difference between g_J, h_J and jac_V
        integer :: id, jd, kd
        
        if(test_this('metric_C2V')) then
            call writo('test whether we have C2V_dn*C2V_up^T = 1')
            call lvl_ud(1)
            
            ! initalize theta and zeta. No need for field-line following theta
            allocate(zeta(n_zeta, n_r)); zeta = 0.0_dp
            allocate(theta(n_zeta, n_r)); theta = 0.0_dp
            do kd = 1,n_r
                zeta(:,kd) = eqd_mesh(n_zeta, min_zeta, max_zeta)
            end do
            theta = zeta
            
            ! calculate the cylindrical variables R and Z and derivatives
            call calc_RZl
            
            ! calculate the metrics in the cylindrical coordinate system
            call metric_C
            
            ! calculate the transformation matrix C(ylindrical) -> V(mec)
            call metric_C2V
            
            ! calculate C2V_dn * C2V_up^T and check if it's equal to 1
            u_mat = 0.0_dp
            u_mat(1,1) = 1.0_dp; u_mat(2,2) = 1.0_dp; u_mat(3,3) = 1.0_dp
            max_index = 0
            diff_max = 0.0_dp
            do kd = 1,n_r
                do jd = 1,n_zeta
                    do id = 1,n_zeta
                        C2V_mult = mat_mult(C2V_up(:,:,id,jd,kd),&
                            &transpose(C2V_dn(:,:,id,jd,kd)))
                        diff_mat = C2V_mult - u_mat
                        if (maxval(abs(diff_mat)).gt.diff_max) then
                            max_index = [id,jd,kd]
                            diff_max = maxval(abs(diff_mat))
                        end if
                    end do
                end do
            end do
            
            call writo('maximum deviation from unity matrix found at ('//&
                &trim(i2str(max_index(1)))//','//trim(i2str(max_index(2)))//&
                &','//trim(i2str(max_index(3)))//'), equal to '//&
                &trim(r2strt(diff_max)))
            
            call lvl_ud(-1)
            
            call writo('Test whether the determinants of g_C, h_C are agreee &
                &with jac_V')
            call lvl_ud(1)
            
            ! calculate  the  metric  factors in the VMEC coordinate system 
            call metric_V
            
            diff_J_max = 0.0_dp
            do kd = 1,n_r
                do jd = 1,n_zeta
                    do id = 1,n_zeta
                        g_J = sqrt(det(3,g_V(:,:,id,jd,kd)))
                        h_J = 1.0_dp/sqrt(det(3,h_V(:,:,id,jd,kd)))
                        if (abs(g_J-jac_V(id,jd,kd)).gt.diff_J_max(1)) then
                            max_J_index(1,:) = [id,jd,kd]
                            diff_J_max(1) = abs(g_J-jac_V(id,jd,kd))
                        end if
                        if (h_J-jac_V(id,jd,kd).gt.diff_J_max(2)) then
                            max_J_index(2,:) = [id,jd,kd]
                            diff_J_max(2) = abs(h_J-jac_V(id,jd,kd)) 
                        end if
                    end do
                end do
            end do
            
            call writo('maximum deviation of sqrt(det(g_J)) from jac_V found &
                & at ('//trim(i2str(max_J_index(1,1)))//','//&
                &trim(i2str(max_J_index(1,2)))//','//&
                &trim(i2str(max_J_index(1,3)))//'), equal to '//&
                &trim(r2strt(diff_J_max(1))))
            call writo('maximum deviation of 1/sqrt(det(h_J)) from jac_V found &
                & at ('//trim(i2str(max_J_index(2,1)))//','//&
                &trim(i2str(max_J_index(2,2)))//','//&
                &trim(i2str(max_J_index(2,3)))//'), equal to '//&
                &trim(r2strt(diff_J_max(2))))
            
            call lvl_ud(-1)
            
            call writo('Paused... press enter')
            read(*,*)
        end if
    end subroutine

    subroutine test_theta_B()
        use VMEC_vars, only: n_r, mpol, ntor,l_c, l_s, iotaf
        use fourier_ops, only: mesh_cs, f2r
        use eq_vars, only: eqd_mesh, theta_B
        use output_ops, only: format_out
        
        integer :: id, kd, n_theta
        real(dp) :: zeta_in(n_r), theta_in(n_r)
        real(dp) :: theta_out(n_r)
        real(dp) :: alpha_in
        real(dp) :: alpha_calc(n_r)
        real(dp) :: cs(0:mpol-1,-ntor:ntor,2)
        real(dp) :: lam(4)
        real(dp), allocatable :: f(:,:), theta_plot(:)
        real(dp), allocatable :: plot_f_theta(:,:)
        integer :: format_out_old
        
        alpha_in = pi*1.2_dp
        zeta_in = pi*0.4_dp
        
        
        if(test_this('theta_B')) then
            call writo('Calculating whether the theta calculated by theta_B &
                &results effectively in the given alpha when substituted')
            call lvl_ud(1)
            
            ! starting value equal to zeta
            theta_in = zeta_in
            theta_out = theta_B(alpha_in,zeta_in,theta_in)
            
            ! calculate lamda to be able to test alpha
            do kd = 1, n_r
                ! cosines and sines
                cs = mesh_cs(mpol,ntor,theta_out(kd),zeta_in(kd))
                
                ! calculate lambda and angular derivatives, converted to FM
                if (kd.eq.1) then                                               ! first point not defined on HM -> extrapolate
                    lam = 3./2. * f2r(l_c(:,:,2),l_s(:,:,2),cs,mpol,ntor) - &
                        &1./2. * f2r(l_c(:,:,3),l_s(:,:,3),cs,mpol,ntor)
                else if (kd.eq.n_r) then                                        ! last point -> extrapolate as well
                    lam = 3./2. * f2r(l_c(:,:,n_r),l_s(:,:,n_r),cs,mpol,ntor)-&
                        &1./2. * f2r(l_c(:,:,n_r-1),l_s(:,:,n_r-1),cs,mpol,ntor)
                else                                                            ! intermediate points -> interpolate
                    lam = 1./2. * f2r(l_c(:,:,kd),l_s(:,:,kd),cs,mpol,ntor) + &
                        &1./2. * f2r(l_c(:,:,kd+1),l_s(:,:,kd+1),cs,mpol,ntor)
                end if
                
                ! alpha
                alpha_calc(kd) = zeta_in(kd) - &
                    &(theta_out(kd) + lam(1))/iotaf(kd)
            end do
            
            call writo('with alpha = '//trim(r2str(alpha_in))//&
                &', and zeta = ')
            call print_ar_1(zeta_in)
            call writo('  -> starting values for theta chosen equal to:')
            call print_ar_1(theta_in)
            call writo('     and final values equal to:')
            call print_ar_1(theta_out)
            call writo('As a check: calculated alpha = ')
            call print_ar_1(alpha_calc)
            call lvl_ud(-1)
            
            call writo('Calculating f for a range of theta values')
            
            call lvl_ud(1)
            n_theta = 50
            
            allocate(theta_plot(n_theta)); theta_plot = 0.0_dp
            theta_plot = eqd_mesh(n_theta, -3_dp*pi, 3_dp*pi)
            allocate(f(n_theta,n_r)); f = 0.0_dp
            allocate(plot_f_theta(2,n_theta))
            
            do kd = 1, n_r
                do id = 1, n_theta
                    ! cosines and sines
                    cs = mesh_cs(mpol,ntor,theta_plot(id),zeta_in(kd))
                    ! lambda
                    lam = f2r(l_c(:,:,kd),l_s(:,:,kd),cs,mpol,ntor)
                    f(id,kd) = zeta_in(kd) - (theta_plot(id)+lam(1))/iotaf(kd) &
                        &- alpha_in
                end do
                plot_f_theta(1,:) = theta_plot
                plot_f_theta(2,:) = f(:,kd)
                
                ! print the output to the screen
                call writo('Plot of f(theta) for r = '//trim(i2str(kd)))
                call print_ar_2(plot_f_theta)
                
                ! plot it on the screen as well using format_out = 3
                format_out_old = format_out
                format_out = 3
                call write_out(2, n_theta, plot_f_theta, &
                    &'f as a function of theta at r = '//trim(i2str(kd)),&
                    &comment='theta_0 = '//trim(r2strt(theta_out(kd)))//&
                    &' for zeta = '//trim(r2strt(zeta_in(kd))))
                format_out = format_out_old
            end do
            
            call writo('Paused... press enter')
            read(*,*)
            
            call lvl_ud(-1)
        end if
    end subroutine

    logical function test_this(test_name)
        use output_ops, only: &
            &lvl_sep
        character(len=*) :: test_name
        character(len=max_str_ln) :: answer_str

        integer :: id 

        call writo('test ' // trim(test_name) // '?')
        do id = 1,lvl                                                           ! to get the response on the right collumn
            write(*,'(A)',advance='no') lvl_sep
        end do
        write(*,'(A)',advance='no') 'y(es)/n(o) [no]: '
        call stop_time
        read (*, '(A)') answer_str
        call start_time

        test_this = .false.
        select case (strh2l(trim(answer_str)))
            !case ('y','Y','yes','Yes','YEs','YES','yEs','yES','yeS','YeS')
            case ('y','yes')
                test_this = .true.
            case ('n','N','no','No','NO','nO') 
            case default 
        end select
    end function test_this
end module test
