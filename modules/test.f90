!-------------------------------------------------------
!   Test some routines and functions
!-------------------------------------------------------
module test
    use num_vars, only: dp, max_str_ln, pi
    use output_ops, only: writo, lvl_ud, print_ar_1, print_ar_2, write_out
    use input_ops, only: yes_no
    use time, only: start_time, stop_time
    use str_ops, only: i2str, r2str, r2strt
    use var_ops, only: mat_mult, det

    implicit none
    private
    public test_repack, test_write_out, test_mesh_cs, test_metric_transf, &
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

    subroutine test_metric_transf()
        use metric_ops, only: metric_C, metric_C2V, metric_V, metric_V2F, &
            &metric_F, &
            &C2V_up, C2V_dn, V2F_up, V2F_dn, jac_V, g_V, h_V, jac_F, g_F, h_F
        use eq_vars, only: eqd_mesh, calc_RZl, calc_flux_q, flux_brkdwn, &
            &n_par, min_par, max_par, theta, zeta
        use VMEC_vars, only: n_r
        
        ! local variables
        integer :: kd
        
        call writo('test metric_transf?')
        if(yes_no(.false.)) then
            ! initalize theta and zeta. No need for field-line following theta
            allocate(zeta(n_par, n_r)); zeta = 0.0_dp
            allocate(theta(n_par, n_r)); theta = 0.0_dp
            do kd = 1,n_r
                zeta(:,kd) = eqd_mesh(n_par, min_par, max_par)
            end do
            theta = zeta
            
            ! calculate the cylindrical variables R, Z and lambda and derivatives
            call calc_RZl
            
            ! calculate flux quantities
            call calc_flux_q
            
            ! find out where using the poloidal flux as normal coordinate breaks
            ! down and calculate the transformation points
            call flux_brkdwn
            
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
            call transf_inv_transf(C2V_up, C2V_dn)
            call lvl_ud(-1)
            call writo('VMEC -> flux:')
            call lvl_ud(1)
            call transf_inv_transf(V2F_up, V2F_dn)
            call lvl_ud(-1)
            
            call lvl_ud(-1)
            
            call writo('Test whether the determinants of g_C, h_C are agreee &
                &with jac_V')
            call lvl_ud(1)
            
            call writo('VMEC:')
            call lvl_ud(1)
            call check_dets(jac_V,g_V,h_V)
            call lvl_ud(-1)
            call writo('flux:')
            call lvl_ud(1)
            call check_dets(jac_F,g_F,h_F)
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
            integer :: max_index(2)                                             ! index of maximum difference
            real(dp) :: diff_max                                                ! maximum deviation of unity matrix
            integer :: id
            
            u_mat = 0.0_dp
            u_mat(1,1) = 1.0_dp; u_mat(2,2) = 1.0_dp; u_mat(3,3) = 1.0_dp
            max_index = 0
            diff_max = 0.0_dp
            do kd = 1,n_r
                do id = 1,n_par
                    T_mult = mat_mult(T_up(:,:,id,kd),&
                        &transpose(T_dn(:,:,id,kd)))
                    diff_mat = T_mult - u_mat
                    if (maxval(abs(diff_mat)).gt.diff_max) then
                        max_index = [id,kd]
                        diff_max = maxval(abs(diff_mat))
                    end if
                end do
            end do
            
            call writo('maximum deviation from unity matrix found at ('//&
                &trim(i2str(max_index(1)))//','//trim(i2str(max_index(2)))//&
                &'), equal to '//trim(r2strt(diff_max)))
        end subroutine
        
        ! check whether jacobians calculated in the routines in metric_ops and
        ! the correct determinants are in agreement
        subroutine check_dets(jac,g,h)
            ! input / output
            real(dp) :: jac(n_par,n_r)
            real(dp) :: g(3,3,n_par,n_r), h(3,3,n_par,n_r)
            
            ! local variables
            real(dp) :: g_J, h_J                                                ! jacobian calculated from g and h
            integer :: max_J_index(2,2)                                         ! index of maximum difference 
            real(dp) :: diff_J_max(2)                                           ! maximum of difference between g_J, h_J and jac_V
            integer :: id, kd
            !real(dp) :: j_calc(n_par,n_r,2)
            
            diff_J_max = 0.0_dp
            do kd = 1,n_r
                do id = 1,n_par
                    g_J = sqrt(abs(det(3,g(:,:,id,kd))))
                    !j_calc(id,kd,1) = g_J
                    h_J = 1.0_dp/sqrt(abs(det(3,h(:,:,id,kd))))
                    !j_calc(id,kd,2) = h_J
                    if (abs(g_J-jac(id,kd)).gt.diff_J_max(1)) then
                        max_J_index(1,:) = [id,kd]
                        diff_J_max(1) = abs(g_J-jac(id,kd))
                    end if
                    if (h_J-jac(id,kd).gt.diff_J_max(2)) then
                        max_J_index(2,:) = [id,kd]
                        diff_J_max(2) = abs(h_J-jac(id,kd))
                    end if
                end do
            end do
            
            call writo('maximum deviation of sqrt(det(g_J)) from jac found &
                & at ('//trim(i2str(max_J_index(1,1)))//','//&
                &trim(i2str(max_J_index(1,2)))//','//&
                &'), equal to '//trim(r2strt(diff_J_max(1))))
            call writo('maximum deviation of 1/sqrt(det(h_J)) from jac found &
                & at ('//trim(i2str(max_J_index(2,1)))//','//&
                &trim(i2str(max_J_index(2,2)))//','//&
                &'), equal to '//trim(r2strt(diff_J_max(2))))
        end subroutine
    end subroutine

    subroutine test_theta_B()
        use VMEC_vars, only: n_r, mpol, ntor
        use eq_vars, only: eqd_mesh, h2f, calc_mesh, &
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
        
        call writo('test theta_B?')
        if(yes_no(.false.)) then
            call writo('Plot zeta(theta)')
            
            call lvl_ud(1)
            
            ! CALCULATE THETA(ZETA)
            call writo('Calculating theta(zeta)')
            alpha = 1.2_dp*pi
            
            ! decrease the number of parallel points for less plots
            n_par_old = n_par
            n_par = 4
            
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
                        plot_var(2,:) = &
                            find_f_plot(n_plot,plot_dep_var,theta(id,kd),alpha) ! corresponding f(zeta)
                    else                                                        ! looking for theta
                        plot_var(1,:) = plot_dep_var                            ! eq. mesh of theta
                        plot_var(2,:) = &
                            find_f_plot(n_plot,plot_dep_var,zeta(id,kd),alpha)  ! corresponding f(theta)
                    end if
                    
                    ! plot it on the screen using format_out = 3
                    format_out_old = format_out
                    format_out = 3
                    call write_out(2,n_plot,plot_var, 'f('//trim(dep_ang)&
                        &//') = zeta -q(theta + lambda) - alpha_0 at r = '//&
                        &trim(i2str(kd))//'/'//trim(i2str(n_r)), comment=&
                        &'= 0 at (theta, zeta) = ('//trim(r2strt(theta(id,kd)))&
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
                    lam = f2r(l_c_F(:,:,kd),l_s_F(:,:,kd),cs,mpol,ntor,1)
                    find_f_plot(jd) = dep_var(jd)-(par_var+lam)/iotaf(kd)-alpha
                else                                                            ! looking for theta (dep. var)
                    ! cosines and sines
                    cs = mesh_cs(mpol,ntor,dep_var(jd),par_var)
                    ! lambda
                    lam = f2r(l_c_F(:,:,kd),l_s_F(:,:,kd),cs,mpol,ntor,1)
                    find_f_plot(jd) = par_var-(dep_var(jd)+lam)/iotaf(kd)-alpha
                end if
            end do
        end function find_f_plot
    end subroutine
end module test
