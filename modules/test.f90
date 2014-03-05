! Test some routines and functions
module test
    use num_vars, only: dp, max_str_ln
    use output_ops, only: writo, lvl_ud
    use var_ops, only: strh2l
    use time, only: start_time, stop_time

    implicit none
    private
    public test_repack, test_write_out, test_mesh_cs, test_calc_metric

contains
    subroutine test_repack()
        ! VMEC variable has structure (1:mnmax, 1:n_r)
        ! output variable should have (1:n_r, 0:mpol-1, -ntor:ntor)
        use fourier_ops, only: repack
        integer :: n_rB, mpolB, ntorB, mnmaxB, nfpB
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
            nfpB = 1
            varinB(:,1) = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
            varinB(:,2) = 2*varinB(:,1)
            varoutB = repack(varinB,mnmaxB,n_rB,mpolB,ntorB,xmB,xnB/nfpB)
                
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
        use plasma_vars, only: &
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
        use plasma_vars, only: mpol, ntor, nfp
     
        real(dp), allocatable :: output(:,:,:)
        real(dp) :: theta, zeta

        allocate(output(0:mpol-1,-ntor:ntor,2))
        theta = 0
        zeta = 0
        
        do
            if(test_this('mesh_cs')) then
                call writo('Testing write_out')
                call lvl_ud(1)
                
                write(*,'(A)',advance='no') 'mpol, ntor, nfp = ' 
                read(*,*) mpol, ntor, nfp
                
                write(*,'(A)',advance='no') 'theta, zeta (*pi) = ' 
                read(*,*) theta, zeta
                theta = pi*theta
                zeta = pi*zeta
                
                output =  mesh_cs(mpol,ntor,nfp,theta,zeta)
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

    subroutine test_calc_metric()
        use plasma_vars, only: calc_metric, n_r, &
            &hrr, hzz, htt, htz, hrt, hrz, grr, gzz, gtt, gtz, grt, grz, jac
        use num_vars, only: n_theta, n_zeta, dp
        use var_ops, only: mat_mult, mat_sub
        use output_ops, only: print_ar_2, write_out
        
        integer :: id, jd, kd
        real(dp) :: g(3,3), h(3,3), g_max(3,3,2), h_max(3,3,2)                  ! lower and upper metric matrix at a 3D point, value at max diff
        real(dp) :: gh(3,3), gh_max(3,3)                                        ! lower metric matrice * upper metric matrix, value at max diff
        real(dp) :: hg(3,3), hg_max(3,3)                                        ! upper metric matrice * lower metric matrix, value at max diff
        real(dp) :: u_mat(3,3), diff_mat(3,3,2)                                 ! unity 3x3 matrix, difference matrix of gh (of hg) - unity_mat
        real(dp) :: diff(n_theta,n_zeta,n_r,2)                                  ! maximum value of difference matrix as function of theta, zeta, r
        integer :: max_index(3,2)                                               ! index of maximum difference
        
        if(test_this('calc_metric')) then
            h = 0.0_dp
            g = 0.0_dp
            u_mat = 0.0_dp
            u_mat(1,1) = 1.0_dp
            u_mat(2,2) = 1.0_dp
            u_mat(3,3) = 1.0_dp
            diff = 0.0_dp
            max_index = 0
            
            do id = 2,n_theta
                do jd = 2,n_zeta
                    do kd = 2, n_r-1                                            ! exclude the values where Jac = 0
                        g(1,1) = grr(id,jd,kd)
                        g(1,2) = grt(id,jd,kd)
                        g(1,3) = grz(id,jd,kd)
                        g(2,1) = grt(id,jd,kd)
                        g(2,2) = gtt(id,jd,kd)
                        g(2,3) = gtz(id,jd,kd)
                        g(3,1) = grz(id,jd,kd)
                        g(3,2) = gtz(id,jd,kd)
                        g(3,3) = gzz(id,jd,kd)
                        h(1,1) = hrr(id,jd,kd)
                        h(1,2) = hrt(id,jd,kd)
                        h(1,3) = hrz(id,jd,kd)
                        h(2,1) = hrt(id,jd,kd)
                        h(2,2) = htt(id,jd,kd)
                        h(2,3) = htz(id,jd,kd)
                        h(3,1) = hrz(id,jd,kd)
                        h(3,2) = htz(id,jd,kd)
                        h(3,3) = hzz(id,jd,kd)
                        gh = mat_mult(g,h)
                        hg = mat_mult(h,g)
                        
                        diff_mat(:,:,1) = abs(mat_sub(u_mat,gh))
                        diff_mat(:,:,2) = abs(mat_sub(u_mat,hg))
                        if (maxval(diff_mat(:,:,1)).gt.maxval(diff(:,:,:,1))) &
                            then
                            max_index(:,1) = [id, jd, kd]
                            g_max(:,:,1) = g
                            h_max(:,:,1) = h
                            gh_max = gh
                        end if
                        diff(id,jd,kd,1) = maxval(diff_mat(:,:,1))
                        if (maxval(diff_mat(:,:,2)).gt.maxval(diff(:,:,:,2))) &
                            then
                            max_index(:,2) = [id, jd, kd]
                            g_max(:,:,2) = g
                            h_max(:,:,2) = h
                            hg_max = hg
                        end if
                        diff(id,jd,kd,2) = maxval(diff_mat(:,:,2))
                    end do
                end do
            end do
            
            write(*,*) 'max diff of gh found at (theta, zeta, r) :', &
                &max_index(:,1)
            write(*,*) 'diff = ', &
                &diff(max_index(1,1),max_index(2,1),max_index(3,1),1)
            write(*,*) 'with g : '
            call print_ar_2(g_max(:,:,1))
            write(*,*) 'and h : '
            call print_ar_2(h_max(:,:,1))
            write(*,*) 'so gh : '
            call print_ar_2(gh_max)
            
            write(*,*) 'At this point, the jacobian is equal to:', &
                &jac(max_index(1,1),max_index(2,1),max_index(3,1))
            
            write(*,*) 'max diff of hg found at (theta, zeta, r) :', &
                &max_index(:,2)
            write(*,*) 'diff = ', &
                &diff(max_index(1,2),max_index(2,2),max_index(3,2),2)
            write(*,*) 'with g : '
            call print_ar_2(g_max(:,:,2))
            write(*,*) 'and h : '
            call print_ar_2(h_max(:,:,2))
            write(*,*) 'so hg : '
            call print_ar_2(hg_max)
            
            write(*,*) 'At this point, the jacobian is equal to:', &
                &jac(max_index(1,2),max_index(2,2),max_index(3,2))
            
            write(*,*) 'differences as function of (theta, zeta, r) have &
                & been written to the output file'
            call write_out(n_theta, n_zeta, diff(:,:,2,1), 'diff_gh_2', &
                &comment = 'diff gh at r=2')
            call write_out(n_theta, n_zeta, diff(:,:,n_r-1,2), 'diff_gh_nr-1', &
                &comment = 'diff gh at r=nr-1')

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
