! Test some routines and functions
module test
    use num_vars, only: dp, max_str_ln
    use output_ops, only: writo, lvl_ud, print_ar_1, print_ar_2
    use var_ops, only: strh2l
    use time, only: start_time, stop_time

    implicit none
    private
    public test_repack, test_write_out, test_mesh_cs, test_metric_C2V

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
        end if
    end subroutine

    subroutine test_write_out()
        use file_ops, only: &
            &output_i
        use output_ops, only: write_out
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
        end if
    end subroutine

    subroutine test_mesh_cs()
        use fourier_ops, only: mesh_cs
        use output_ops, only: print_ar_2
        use num_vars, only: pi
        use VMEC_vars, only: mpol, ntor
     
        real(dp), allocatable :: output(:,:,:)
        real(dp) :: theta, zeta

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
        end do
    end subroutine

    subroutine test_metric_C2V()
        use num_vars, only: pi
        use metric_ops, only: metric_C, metric_C2V, metric_V, &
            &C2V_up, C2V_dn, jac_V, g_V, h_V
        use grid_vars, only: calc_ang_mesh, calc_RZl, &
            &n_theta, n_zeta, theta, zeta
        use VMEC_vars, only: n_r
        use var_ops, only: i2str, r2strt, mat_mult, det
        
        real(dp) :: min_theta, max_theta, min_zeta, max_zeta
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
            
            n_theta = 10; min_theta = 0; max_theta = 2*pi
            n_zeta = 10; min_zeta = 0; max_zeta = 2*pi
            theta = calc_ang_mesh(n_theta, min_theta, max_theta)
            zeta = calc_ang_mesh(n_zeta, min_zeta, max_zeta)
            
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
                    do id = 1,n_theta
                        C2V_mult = mat_mult(C2V_up(id,jd,kd,:,:),&
                            &transpose(C2V_dn(id,jd,kd,:,:)))
                        diff_mat = C2V_mult - u_mat
                        if (maxval(diff_mat).gt.diff_max) then
                            max_index = [id,jd,kd]
                            diff_max = maxval(diff_mat)
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
                    do id = 1,n_theta
                        g_J = sqrt(det(3,g_V(id,jd,kd,:,:)))
                        h_J = 1.0_dp/sqrt(det(3,h_V(id,jd,kd,:,:)))
                        if (g_J-jac_V(id,jd,kd).gt.diff_J_max(1)) then
                            max_J_index(1,:) = [id,jd,kd]
                            diff_J_max(1) = g_J-jac_V(id,jd,kd) 
                        end if
                        if (h_J-jac_V(id,jd,kd).gt.diff_J_max(2)) then
                            max_J_index(2,:) = [id,jd,kd]
                            diff_J_max(2) = h_J-jac_V(id,jd,kd) 
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
            
        end if
    end subroutine

    logical function test_this(test_name)
        use output_ops, only: lvl, lvl_sep
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
