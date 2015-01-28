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
    public test_repack, test_print_GP, test_metric_transf, test_ang_B, &
        &test_calc_ext_var, test_B, test_calc_deriv, test_add_arr_mult, &
        &test_conv_FHM, test_calc_RZL, test_calc_T_VF, test_calc_inv_met, &
        &test_calc_det, test_inv, test_calc_f_deriv, test_calc_g, &
        &test_fourier2real, test_prepare_X, test_slepc, test_pres_balance
    
contains
    integer function test_fourier2real() result(ierr)
        !use VMEC, only: nfp
        !use eq_vars, only: calc_eqd_grid
        !use fourier_ops, only: fourier2real, calc_trigon_factors
        
        !character(*), parameter :: rout_name = 'test_fourier2real'
        
        !! local variables
        !real(dp), allocatable :: cosvar(:,:)
        !real(dp), allocatable :: sinvar(:,:)
        !real(dp), allocatable :: realvar(:,:)
        !real(dp), allocatable :: realvar_ALT(:)
        !real(dp), allocatable :: theta(:,:), zeta(:,:)
        !real(dp), allocatable :: trigon_factors(:,:,:,:,:)                      ! trigonometric factor cosine for the inverse fourier transf.
        !integer :: id, jd
        !real(dp) :: min_n, max_n
        !integer :: ntor, mpol
        !integer :: n
        
        !! initialize ierr
        !ierr = 0
        
        !write(*,*) 'test fourier2real?'
        !if(yes_no(.false.)) then
            
            !output_i = 0
            !mpol = 3
            !ntor = 1
            !n = 200
            
            !allocate(cosvar(0:mpol-1,-ntor:ntor))
            !allocate(sinvar(0:mpol-1,-ntor:ntor))
            !allocate(realvar(n,1))
            !allocate(realvar_ALT(n))
            
            !allocate(theta(n,1))
            !allocate(zeta(n,1))
            
            !min_n = 0.0_dp
            !max_n = 3*pi
            
            !! physical grid
            !zeta = 0.4*pi/2
            !ierr = calc_eqd_grid(theta(:,1),n,min_n,max_n)
            !CHCKERR('')
            
            !! set up Fourier coefficients
            !cosvar(:,-1) = [3,2,3]
            !sinvar(:,-1) = [0,1,0]
            !cosvar(:,0) = [4,5,0]
            !sinvar(:,0) = [1,0,3]
            !cosvar(:,1) = [1,2,1]
            !sinvar(:,1) = [1,3,4]
            
            !write(*,*) 'zeroth order'
            !! invert Fourier series
            !ierr = calc_trigon_factors(theta,zeta,trigon_factors)
            !CHCKERR('')
            !ierr = fourier2real(cosvar,sinvar,trigon_factors,&
                !&realvar)
            
            !call print_GP_2D('real variable','',realvar)
            
            !! do an alternative calculation
            !realvar_ALT = 0.0_dp
            !do id = 0,mpol-1
                !do jd = -ntor,ntor
                    !realvar_ALT = realvar_ALT + &
                        !&cosvar(id,jd)*cos(id*theta-jd*nfp*zeta) + &
                        !&sinvar(id,jd)*sin(id*theta-jd*nfp*zeta)
                !end do
            !end do
            
            !call print_GP_2D('real variable, alternative calc','',realvar_ALT)
            
            !! difference
            !call print_GP_2D('difference','',realvar-realvar_ALT)
            
            !write(*,*) 'D_theta'
            !! invert Fourier series
            !do id = 1,n
                !ierr = calc_ang_grid(cs,mpol,ntor,nfp,theta(id),zeta(id))
                !CHCKERR('')
                !realvar(id) = f2r(cosvar,sinvar,cs,mpol,ntor,nfp,[1,0],ierr)
                !CHCKERR('')
            !end do
            
            !call print_GP_2D('real variable','',realvar)
            
            !! do an alternative calculation
            !realvar_ALT = 0.0_dp
            !do id = 0,mpol-1
                !do jd = -ntor,ntor
                    !realvar_ALT = realvar_ALT + &
                        !&id*sinvar(id,jd)*cos(id*theta-jd*nfp*zeta) - &
                        !&id*cosvar(id,jd)*sin(id*theta-jd*nfp*zeta)
                !end do
            !end do
            
            !call print_GP_2D('real variable, alternative calc.','',realvar_ALT)
            
            !! difference
            !call print_GP_2D('difference','',realvar-realvar_ALT)
            
            !write(*,*) 'D_zeta'
            !! invert Fourier series
            !do id = 1,n
                !ierr = calc_ang_grid(cs,mpol,ntor,nfp,theta(id),zeta(id))
                !CHCKERR('')
                !realvar(id) = f2r(cosvar,sinvar,cs,mpol,ntor,nfp,[0,1],ierr)
                !CHCKERR('')
            !end do
            
            !call print_GP_2D('real variable','',realvar)
            
            !! do an alternative calculation
            !realvar_ALT = 0.0_dp
            !do id = 0,mpol-1
                !do jd = -ntor,ntor
                    !realvar_ALT = realvar_ALT - &
                        !&jd*sinvar(id,jd)*cos(id*theta-jd*nfp*zeta) + &
                        !&jd*cosvar(id,jd)*sin(id*theta-jd*nfp*zeta)
                !end do
            !end do
            
            !call print_GP_2D('real variable, alternative calc.','',realvar_ALT)
            
            !! difference
            !call print_GP_2D('difference','',realvar-realvar_ALT)
            
            !write(*,*) 'D^2_theta D_zeta'
            !! invert Fourier series
            !do id = 1,n
                !ierr = calc_ang_grid(cs,mpol,ntor,nfp,theta(id),zeta(id))
                !CHCKERR('')
                !realvar(id) = f2r(cosvar,sinvar,cs,mpol,ntor,nfp,[2,1],ierr)
                !CHCKERR('')
            !end do
            
            !call print_GP_2D('real variable','',realvar)
            
            !! do an alternative calculation
            !realvar_ALT = 0.0_dp
            !do id = 0,mpol-1
                !do jd = -ntor,ntor
                    !realvar_ALT = realvar_ALT + &
                        !&id*id*jd*sinvar(id,jd)*cos(id*theta-jd*nfp*zeta) - &
                        !&id*id*jd*cosvar(id,jd)*sin(id*theta-jd*nfp*zeta)
                !end do
            !end do
            
            !call print_GP_2D('real variable, alternative calc.','',realvar_ALT)
            
            !! difference
            !call print_GP_2D('difference','',realvar-realvar_ALT)
            
            
            !write(*,*) 'Stopping'
            !stop
        !end if
    end function test_fourier2real
    
    subroutine test_repack
        ! VMEC variable has structure (1:mnmax, 1:grp_n_r_eq)
        ! output variable should have (1:grp_n_r_eq, 0:mpol-1, -ntor:ntor)
        use VMEC, only: repack
        
        integer :: n_rB, mpolB, ntorB, mnmaxB
        real(dp), allocatable :: xmB(:), xnB(:)
        real(dp), allocatable :: varinB(:,:)
        real(dp), allocatable :: varoutB(:,:,:)
            
        write(*,*) 'test repack?'
        if(yes_no(.false.)) then
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
            
            write(*,*) 'Paused... press enter'
            read(*,*)
            
            write(*,*) 'Stopping'
            stop
        end if
    end subroutine test_repack

    subroutine test_print_GP
        use output_ops, only: print_GP_2D, print_GP_3D
        
        ! local variables
        integer :: test_type
        integer :: nplt, npnt, npntx, npnty
        real(dp), allocatable :: y_1(:)
        real(dp), allocatable :: y_2(:,:)
        real(dp), allocatable :: y_3(:,:,:)
        real(dp), allocatable :: x_1(:)
        real(dp), allocatable :: x_2(:,:)
        real(dp), allocatable :: x_3(:,:,:)
        real(dp), allocatable :: z_2(:,:)
        real(dp), allocatable :: z_3(:,:,:)
        integer :: iplt, ipnt, ipntx, ipnty
        
        write(*,*) 'test print_GP?'
        if(yes_no(.false.)) then
            write(*,*) 'Testing print_GP'
            
            do
                write(*,*) 'Which type of test do you want to perform?'
                write(*,*) '    2: 2D array'
                write(*,*) '    3: 3D array'
                read(*,*) test_type
                if (test_type.ge.2 .and. test_type.le.3) then
                    exit
                else
                    write(*,*) 'please select a valid test type'
                end if
            end do
            
            select case (test_type)
                case(2)
                    call test_2D
                case(3)
                    call test_3D
                case default
                    write(*,*) 'ERROR: In test_print_GP, no test type &
                        &associated with '//trim(i2str(test_type))//'...'
                    stop
            end select
            
            
            
            write(*,*) 'Stopping'
            stop
        end if
    contains
        subroutine test_2D
            ! one function, no x values
            write(*,*) 'plotting data without x values'
            ! allocate
            npnt = 20
            allocate(y_1(npnt))
            ! set up some data
            do ipnt = 1,npnt
                y_1(ipnt) = cos((ipnt-1)*3*pi/(npnt-1))
            end do
            ! plot this data
            call print_GP_2D('cos(3*x)','test1',y_1)
            ! deallocate
            deallocate(y_1)
            
            ! one function, with x values
            write(*,*) 'plotting data with x values'
            ! allocate
            npnt = 20
            allocate(x_1(npnt))
            allocate(y_1(npnt))
            ! set up some data
            do ipnt = 1,npnt
                x_1(ipnt) = 5*ipnt
                y_1(ipnt) = cos((ipnt-1)*3*pi/(npnt-1))
            end do
            ! plot this data
            call print_GP_2D('cos(3*x) with x','test2',y_1,x_1)
            ! deallocate
            deallocate(x_1)
            deallocate(y_1)
            
            ! many functions, without x values
            write(*,*) 'plotting multiple data without x values'
            ! allocate
            npnt = 20
            nplt = 5
            allocate(y_2(npnt,nplt))
            ! set up some data
            do iplt = 1,nplt
                do ipnt = 1,npnt
                    y_2(ipnt,iplt) = cos(1.0_dp*iplt/nplt*(ipnt-1)*3*pi/&
                        &(npnt-1))
                end do
            end do
            ! plot this data
            call print_GP_2D('cos(i*3*x)','',y_2)
            ! deallocate
            deallocate(y_2)
            
            ! many functions, with x values
            write(*,*) 'plotting multiple data with x values'
            ! allocate
            npnt = 20
            nplt = 5
            allocate(x_2(npnt,nplt))
            allocate(y_2(npnt,nplt))
            ! set up some data
            do iplt = 1,nplt
                do ipnt = 1,npnt
                    x_2(ipnt,iplt) = 5*ipnt
                    y_2(ipnt,iplt) = cos(1.0_dp*iplt/nplt*(ipnt-1)*3*pi/&
                        &(npnt-1))
                end do
            end do
            ! plot this data
            call print_GP_2D('cos(i*3*x) with x','',y_2,x_2)
            ! deallocate
            deallocate(x_2)
            deallocate(y_2)
        end subroutine
            
        subroutine test_3D
            ! one function, no x values
            write(*,*) 'plotting data without x and y values'
            ! allocate
            npntx = 20
            npnty = 30
            allocate(z_2(npntx,npnty))
            ! set up some data
            do ipntx = 1,npntx
                do ipnty = 1, npnty
                    z_2(ipntx,ipnty) = sin((ipntx-1)*pi/(npntx-1))*&
                        &cos((ipnty-1)*3*pi/(npnty-1))
                end do
            end do
            ! plot this data
            call print_GP_3D('sin(x)*cos(3*y)','',z_2)
            ! deallocate
            deallocate(z_2)
            
            ! one function, with x values
            write(*,*) 'plotting data with x, without y values'
            ! allocate
            npntx = 20
            npnty = 30
            allocate(x_2(npntx,npnty))
            allocate(z_2(npntx,npnty))
            ! set up some data
            do ipntx = 1,npntx
                x_2(ipntx,:) = 5*ipntx
                do ipnty = 1, npnty
                    z_2(ipntx,ipnty) = sin((ipntx-1)*pi/(npntx-1))*&
                        &cos((ipnty-1)*3*pi/(npnty-1))
                end do
            end do
            ! plot this data
            call print_GP_3D('sin(x)*cos(3*y) with x','',z_2,x=x_2)
            ! deallocate
            deallocate(x_2)
            deallocate(z_2)
            
            ! one function, with y values
            write(*,*) 'plotting data with y, without x values'
            ! allocate
            npntx = 20
            npnty = 30
            allocate(y_2(npntx,npnty))
            allocate(z_2(npntx,npnty))
            ! set up some data
            do ipnty = 1,npnty
                y_2(:,ipnty) = 5*ipnty
                do ipntx = 1, npntx
                    z_2(ipntx,ipnty) = sin((ipntx-1)*pi/(npntx-1))*&
                        &cos((ipnty-1)*3*pi/(npnty-1))
                end do
            end do
            ! plot this data
            call print_GP_3D('sin(x)*cos(3*y) with y','',z_2,y=y_2)
            ! deallocate
            deallocate(y_2)
            deallocate(z_2)
            
            ! many functions, without x and y values
            write(*,*) 'plotting multiple data without x and y values'
            ! allocate
            nplt = 3
            npntx = 20
            npnty = 30
            allocate(z_3(npntx,npnty,nplt))
            ! set up some data
            do iplt = 1,nplt
                do ipntx = 1,npntx
                    do ipnty = 1,npnty
                        z_3(ipntx,ipnty,iplt) = sin((ipntx-1)*pi/(npntx-1))*&
                            &cos((ipnty-1)*3*pi/(npnty-1))*1.0_dp*iplt/nplt
                    end do
                end do
            end do
            ! plot this data
            call print_GP_3D('sin(x)*cos(3*y)*i','',z_3)
            ! deallocate
            deallocate(z_3)
            
            !! many functions, with x values
            write(*,*) 'plotting multiple data with x values'
            ! allocate
            nplt = 3
            npntx = 20
            npnty = 30
            allocate(z_3(npntx,npnty,nplt))
            allocate(x_3(npntx,npnty,nplt))
            ! set up some data
            do iplt = 1,nplt
                do ipntx = 1,npntx
                    x_3(ipntx,:,iplt) = 5*ipntx
                    do ipnty = 1,npnty
                        z_3(ipntx,ipnty,iplt) = sin((ipntx-1)*pi/(npntx-1))*&
                            &cos((ipnty-1)*3*pi/(npnty-1))*1.0_dp*iplt/nplt
                    end do
                end do
            end do
            ! plot this data
            call print_GP_3D('sin(x)*cos(3*y)*i with x','',z_3,x=x_3)
            ! deallocate
            deallocate(z_3)
            deallocate(x_3)
            
            !! many functions, with y values
            write(*,*) 'plotting multiple data with y values'
            ! allocate
            nplt = 3
            npntx = 20
            npnty = 30
            allocate(z_3(npntx,npnty,nplt))
            allocate(y_3(npntx,npnty,nplt))
            ! set up some data
            do iplt = 1,nplt
                do ipnty = 1,npnty
                    y_3(:,ipnty,iplt) = 5*ipnty
                    do ipntx = 1,npntx
                        z_3(ipntx,ipnty,iplt) = sin((ipntx-1)*pi/(npntx-1))*&
                            &cos((ipnty-1)*3*pi/(npnty-1))*1.0_dp*iplt/nplt
                    end do
                end do
            end do
            ! plot this data
            call print_GP_3D('sin(x)*cos(3*y)*i with y','random_test_name.txt',&
                &z_3,y=y_3)
            ! deallocate
            deallocate(z_3)
            deallocate(y_3)
            
            !! many functions, with x and y values
            write(*,*) 'plotting multiple data with x and y values'
            ! allocate
            nplt = 3
            npntx = 20
            npnty = 30
            allocate(z_3(npntx,npnty,nplt))
            allocate(y_3(npntx,npnty,nplt))
            allocate(x_3(npntx,npnty,nplt))
            ! set up some data
            do iplt = 1,nplt
                do ipnty = 1,npnty
                    y_3(:,ipnty,iplt) = 5*ipnty
                    do ipntx = 1,npntx
                        z_3(ipntx,ipnty,iplt) = sin((ipntx-1)*pi/(npntx-1))*&
                            &cos((ipnty-1)*3*pi/(npnty-1))*1.0_dp*iplt/nplt
                    end do
                end do
                do ipntx = 1,npntx
                    x_3(ipntx,:,iplt) = 5*ipntx
                end do
            end do
            ! plot this data
            call print_GP_3D('sin(x)*cos(3*y)*i with x and y','',z_3,x=x_3,&
                &y=y_3)
            ! deallocate
            deallocate(z_3)
            deallocate(y_3)
            deallocate(x_3)
        end subroutine
    end subroutine test_print_GP

    integer function test_prepare_X() result(ierr)
        use X_ops, only: prepare_X
        use X_vars, only: PV0, PV1, PV2
        use driver_rich, only: calc_eq
        
        character(*), parameter :: rout_name = 'test_prepare_X'
        
        ! initialize ierr
        ierr = 0
        
        write(*,*) 'test prepare_X?'
        if(yes_no(.false.)) then
            ! calculate equilibrium
            write(*,*) 'calculating equilibrium'
            ierr = calc_eq(0.2*pi)
            CHCKERR('')
            
            ! calculate equilibrium
            write(*,*) 'calculating P'
            ierr = prepare_X()
            
            ! visualize PV
            write(*,*) 'Real PV0 (5,5) ='
            call print_ar_2(realpart(PV0(5,5,:,:)))
            write(*,*) 'Imag PV0 (5,5) ='
            call print_ar_2(imagpart(PV0(5,5,:,:)))
            
            write(*,*) 'Real PV1 (5,5) ='
            call print_ar_2(realpart(PV1(5,5,:,:)))
            write(*,*) 'Imag PV1 (5,5) ='
            call print_ar_2(imagpart(PV1(5,5,:,:)))
            
            write(*,*) 'Real PV2 (5,5) ='
            call print_ar_2(realpart(PV2(5,5,:,:)))
            write(*,*) 'Imag PV2 (5,5) ='
            call print_ar_2(imagpart(PV2(5,5,:,:)))
            
            write(*,*) 'Stopping'
            stop
        end if
    end function test_prepare_X
    
    integer function test_calc_deriv() result(ierr)
        use utilities, only: calc_deriv
        
        character(*), parameter :: rout_name = 'test_calc_deriv'
        
        ! local variables
        integer :: loc_max, id, kd
        integer :: n_steps, grid_type
        real(dp) :: x_loc
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
    
    integer function test_calc_ext_var() result(ierr)
        use utilities, only : calc_ext_var
        
        character(*), parameter :: rout_name = 'test_calc_ext_var'
        
        ! local variables
        integer :: jd, kd
        integer :: n_r
        real(dp), allocatable :: x(:), varin(:), var_num(:),varder(:), &
            &varder_num(:), x_int(:)
        
        ! initialize ierr
        ierr = 0
        
        write(*,*) 'test calc_ext_var'
        if(yes_no(.false.)) then
            output_i = 0
            write(*,*) 'Testing whether extrapolation is working correctly'
            
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
                ierr = calc_ext_var(var_num(kd),[(varin(jd),jd=kd,kd+2)],&
                    &[(x(jd),jd=kd,kd+2)],x_int(kd),0)
                CHCKERR('')
                ierr = calc_ext_var(varder_num(kd),[(varin(jd),jd=kd,kd+2)],&
                    &[(x(jd),jd=kd,kd+2)],x(kd),1)
                CHCKERR('')
            end do
            write(*,*) 'Don''t worry about the last 3 points!!!'
            call print_GP_2D('function: sin(x)+0.5*cos(3x)','',varin,x=x)
            call print_GP_2D('num interp. of function: sin(x)+0.5*cos(3x)','',&
                &var_num,x=x_int)
            call print_GP_2D('deriv of function: cos(x)-1.5*sin(3x)','',&
                &varder,x=x)
            call print_GP_2D('num interp. of deriv of function: &
                &cos(x)-1.5*sin(3x)','',varder_num,x=x)
            call print_GP_2D('relative error: cos(x)-1.5*sin(3x)','',&
                &abs(varder_num-varder),x=x)
            
            !! check whether the interpolating polynomial yields the source points
            !do kd = 1,n_r-3
                !ierr = calc_ext_var(var_num(kd),[(varin(jd),jd=kd,kd+2)],&
                    !&[(x(jd),jd=kd,kd+2)],x(kd),0)
                CHCKERR('')
                !write(*,*) 'kd = ', kd
                !write(*,*) 'var_num = ', var_num(kd)
                !write(*,*) 'varin   = ', varin(kd)
            !end do
            
            write(*,*) 'Stopping'
            stop
        end if
    end function test_calc_ext_var
    
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
    
    integer function test_calc_RZL() result(ierr)
        use coord_ops, only: calc_eqd_grid
        use eq_ops, only: calc_RZL, init_eq
        use eq_vars, only: R_E, Z_E, theta_E, zeta_E,  grp_n_r_eq
        use X_vars, only: n_par
        use VMEC, only: rmax_surf, rmin_surf, zmax_surf
        
        character(*), parameter :: rout_name = 'test_calc_RZL'
        
        ! local variables
        integer :: jd
        
        ! initialize ierr
        ierr = 0
        
        write(*,*) 'test calc_RZL?'
        if(yes_no(.false.)) then
            output_i = 0
            write(*,*) 'initializing'
            ierr = init_eq()
            
            write(*,*) 'calculations with constant zeta, varying theta_E'
            if (allocated(zeta_E)) deallocate(zeta_E)
            allocate(zeta_E(n_par,grp_n_r_eq)); zeta_E = 0.0_dp
            if (allocated(theta_E)) deallocate(theta_E)
            allocate(theta_E(n_par,grp_n_r_eq)); theta_E = 0.0_dp
            
            zeta_E = 0.4*pi/2
            do jd = 1,grp_n_r_eq
                ierr = calc_eqd_grid(theta_E(:,jd),n_par, 0.0_dp*pi, 3.0_dp*pi)
                CHCKERR('')
            end do
            
            ierr = calc_RZL([0,0,0])
            CHCKERR('')
            ierr = calc_RZL([0,1,0])
            CHCKERR('')
            ierr = calc_RZL([0,2,0])
            CHCKERR('')
            ierr = calc_RZL([0,3,0])
            CHCKERR('')
            ierr = calc_RZL([0,0,1])
            CHCKERR('')
            ierr = calc_RZL([0,0,2])
            CHCKERR('')
            ierr = calc_RZL([0,0,3])
            CHCKERR('')
            
            write(*,*) 'Plotting R'
            call print_GP_3D('R','',R_E(:,:,0,0,0))
            call print_GP_2D('R at r = 5','',R_E(:,5,0,0,0),x=theta_E(:,5))
            do jd = 1,3
                call print_GP_2D('R_theta_E^'//trim(i2str(jd))//' at r = 5','',&
                    &R_E(:,5,0,jd,0),x=theta_E(:,5))
            end do
            do jd = 1,3
                call print_GP_2D('R_zeta^'//trim(i2str(jd))//' at r = 5','',&
                    &R_E(:,5,0,0,jd),x=theta_E(:,5))
            end do
            
            write(*,*) 'Plotting Z'
            call print_GP_3D('Z','',Z_E(:,:,0,0,0))
            call print_GP_2D('Z at r = 5','',Z_E(:,5,0,0,0),x=theta_E(:,5))
            do jd = 1,3
                call print_GP_2D('Z_theta_E^'//trim(i2str(jd))//' at r = 5','',&
                    &Z_E(:,5,0,jd,0),x=theta_E(:,5))
            end do
            do jd = 1,3
                call print_GP_2D('Z_zeta^'//trim(i2str(jd))//' at r = 5','',&
                    &Z_E(:,5,0,0,jd),x=theta_E(:,5))
            end do
            
            ! test whether VMEC given bounds are respected
            if (rmax_surf.gt.0.0_dp) then
                write(*,*) 'Checking the bounds of R'
                call within_bounds(R_E(:,:,0,0,0),rmin_surf,rmax_surf)
            else
                write(*,*) 'Not possible to check the bounds of R'
            end if
            if (zmax_surf.gt.0.0_dp) then
                write(*,*) 'Checking the bounds of Z'
                call within_bounds(Z_E(:,:,0,0,0),-zmax_surf,zmax_surf)
            else
                write(*,*) 'Not possible to check the bounds of Z'
            end if
            
            write(*,*) 'Stopping'
            stop
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
                write(*,*) 'WARNING: minimum of variable in real angular space &
                    & is lower than VMEC provided minimum by '//&
                    &trim(r2strt(100*min_frac))//'%...'
                write(*,*) ' -> Maybe use an even number of poloidal points, &
                    &or more angular grid points in general?'
                write(*,*) ' -> Maybe run VMEC with more accuracy?'
            else if (max_frac.gt.margin) then                                   ! too high maximum
                write(*,*) 'WARNING: maximum of variable in real angular space &
                    & is greater than VMEC provided maximum by '//&
                    &trim(r2strt(100*max_frac))//'%...'
                write(*,*) ' -> Maybe use more angular grid points'
                write(*,*) ' -> Maybe run VMEC with more accuracy?'
            else 
                write(*,*) 'within bounds'
                return
            end if
        end subroutine
    end function test_calc_RZL
    
    integer function test_calc_T_VF() result(ierr)
        use driver_rich, only: calc_eq
        use metric_vars, only: T_EF
        use eq_vars, only: q_saf_E, L_E, flux_p_E, theta => theta_E, &
            &grp_n_r_eq
        
        character(*), parameter :: rout_name = 'test_calc_T_VF'
        
        ! local variables
        real(dp) :: T_loc(1:grp_n_r_eq,3,3)
        integer :: id, jd
        integer :: case_nr
        integer :: deriv(3)
        
        ! initialize ierr
        ierr = 0
        
        write(*,*) 'test calc_T_VF?'
        if(yes_no(.false.)) then
            output_i = 0
            do 
                write(*,*) 'Which case do you want?'
                write(*,*) '    1: 0th order'
                write(*,*) '    2: 1st order, r derivative'
                write(*,*) '    3: 1st order, theta derivative'
                write(*,*) '    4: 1st order, zeta derivative'
                write(*,*) '    5: 4th order, rthetathetazeta derivative'
                read(*,*) case_nr
                if(case_nr.lt.1 .or. case_nr.gt.5) then
                    write(*,*) 'You have to choose from the available options'
                else
                    exit
                end if
            end do
            
            write(*,*) 'calculating equilibrium'
            ierr = calc_eq(0.2*pi)
            CHCKERR('')
            
            write(*,*) 'calculating case'
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
                    write(*,*) 'ERROR: In test_T_VF, case '//&
                        &trim(i2str(case_nr))//' is not associated with &
                        &anything. How did you get here???'
                    stop
            end select
            write(*,*) 'case calculated'
            
            write(*,*) 'plotting case'
            do id = 1,3
                do jd = 1,3
                    call print_GP_2D('T_VF('//trim(i2str(id))//','//&
                        &trim(i2str(jd))//')','',T_EF(5,:,id,jd,deriv(1),&
                        &deriv(2),deriv(3)))
                    call print_GP_2D('alt. T_VF('//trim(i2str(id))//','//&
                        &trim(i2str(jd))//')','',T_loc(:,id,jd))
                    call print_GP_2D('abs. diff ('//trim(i2str(id))//','//&
                        &trim(i2str(jd))//')','',abs(T_loc(:,id,jd)-&
                        &T_EF(5,:,id,jd,deriv(1),&
                        &deriv(2),deriv(3))))
                end do
            end do
            write(*,*) 'case plotted'
            
            write(*,*) 'Stopping'
            stop
        end if
    contains
        subroutine case_1
            T_loc = 0.0_dp
            T_loc(:,1,1) = -q_saf_E(:,1)*(theta(5,:)+L_E(5,:,0,0,0)) &
                &-q_saf_E(:,0)*L_E(5,:,1,0,0)
            T_loc(:,1,2) = flux_p_E(:,1)/(2*pi)
            T_loc(:,1,3) = L_E(5,:,1,0,0)
            T_loc(:,2,1) = -q_saf_E(:,0)*(1.0_dp+L_E(5,:,0,1,0))
            T_loc(:,2,2) = 0.0_dp
            T_loc(:,2,3) = 1.0_dp + L_E(5,:,0,1,0)
            T_loc(:,3,1) = 1.0_dp - q_saf_E(:,0)*L_E(5,:,0,0,1)
            T_loc(:,3,2) = 0.0_dp
            T_loc(:,3,3) = L_E(5,:,0,0,1)
        end subroutine
        
        subroutine case_2
            T_loc = 0.0_dp
            T_loc(:,1,1) = -q_saf_E(:,2)*(theta(5,:)+L_E(5,:,0,0,0)) &
                &-2*q_saf_E(:,1)*L_E(5,:,1,0,0)&
                &-q_saf_E(:,0)*L_E(5,:,2,0,0)
            T_loc(:,1,2) = flux_p_E(:,2)/(2*pi)
            T_loc(:,1,3) = L_E(5,:,2,0,0)
            T_loc(:,2,1) = -q_saf_E(:,1)*(1.0_dp+L_E(5,:,0,1,0)) &
                &-q_saf_E(:,0)*L_E(5,:,1,1,0)
            T_loc(:,2,2) = 0.0_dp
            T_loc(:,2,3) = L_E(5,:,1,1,0)
            T_loc(:,3,1) = -q_saf_E(:,1)*L_E(5,:,0,0,1) &
                &- q_saf_E(:,0)*L_E(5,:,1,0,1)
            T_loc(:,3,2) = 0.0_dp
            T_loc(:,3,3) = L_E(5,:,1,0,1)
        end subroutine
        
        subroutine case_3
            T_loc = 0.0_dp
            T_loc(:,1,1) = -q_saf_E(:,1)*(1+L_E(5,:,0,1,0))&
                &-q_saf_E(:,0)*L_E(5,:,1,1,0)
            T_loc(:,1,2) = 0.0_dp
            T_loc(:,1,3) = L_E(5,:,1,1,0)
            T_loc(:,2,1) = -q_saf_E(:,0)*L_E(5,:,0,2,0)
            T_loc(:,2,2) = 0.0_dp
            T_loc(:,2,3) = L_E(5,:,0,2,0)
            T_loc(:,3,1) = -q_saf_E(:,0)*L_E(5,:,0,1,1)
            T_loc(:,3,2) = 0.0_dp
            T_loc(:,3,3) = L_E(5,:,0,1,1)
        end subroutine
        
        subroutine case_4
            T_loc = 0.0_dp
            T_loc(:,1,1) = -q_saf_E(:,1)*L_E(5,:,0,0,1)&
                &-q_saf_E(:,0)*L_E(5,:,1,0,1)
            T_loc(:,1,2) = 0.0_dp
            T_loc(:,1,3) = L_E(5,:,1,0,1)
            T_loc(:,2,1) = -q_saf_E(:,0)*L_E(5,:,0,1,1)
            T_loc(:,2,2) = 0.0_dp
            T_loc(:,2,3) = L_E(5,:,0,1,1)
            T_loc(:,3,1) = -q_saf_E(:,0)*L_E(5,:,0,0,2)
            T_loc(:,3,2) = 0.0_dp
            T_loc(:,3,3) = L_E(5,:,0,0,2)
        end subroutine
        
        subroutine case_5
            T_loc = 0.0_dp
            T_loc(:,1,1) = -q_saf_E(:,2)*L_E(5,:,0,2,1)-2*q_saf_E(:,1)*&
                &L_E(5,:,1,2,1)-q_saf_E(:,0)*L_E(5,:,2,2,1)
            T_loc(:,1,2) = 0.0_dp
            T_loc(:,1,3) = L_E(5,:,2,2,1)
            T_loc(:,2,1) = -q_saf_E(:,1)*L_E(5,:,0,3,1)&
                &-q_saf_E(:,0)*L_E(5,:,1,3,1)
            T_loc(:,2,2) = 0.0_dp
            T_loc(:,2,3) = L_E(5,:,1,3,1)
            T_loc(:,3,1) = -q_saf_E(:,1)*L_E(5,:,0,2,2)&
                &-q_saf_E(:,0)*L_E(5,:,1,2,2)
            T_loc(:,3,2) = 0.0_dp
            T_loc(:,3,3) = L_E(5,:,1,2,2)
        end subroutine
    end function test_calc_T_VF
    
    integer function test_calc_g() result(ierr)
        use driver_rich, only: calc_eq
        use metric_vars, only: g_E, g_F
        use eq_vars, only: R_E, Z_E, q_saf_E
        
        character(*), parameter :: rout_name = 'tset_calc_g'
        
        ! local variables
        real(dp) :: g(3,3)
        integer :: id_d
        integer :: d(3)
        
        ! initialize ierr
        ierr = 0
        
        write(*,*) 'test calc_g?'
        if(yes_no(.false.)) then
            
            write(*,*) 'calculating equilibrium'
            ierr = calc_eq(0.2*pi)
            CHCKERR('')
            
            write(*,*) 'check g_V?'
            if(yes_no(.false.)) then
                write(*,*) 'g_V(5,5) = '
                call print_ar_2(g_E(5,5,:,:,0,0,0))
                
                write(*,*) 'g_V(5,5), manually'
                g(1,1) = R_E(5,5,1,0,0)*R_E(5,5,1,0,0) + &
                    &Z_E(5,5,1,0,0)*Z_E(5,5,1,0,0)
                g(1,2) = R_E(5,5,1,0,0)*R_E(5,5,0,1,0) + &
                    &Z_E(5,5,1,0,0)*Z_E(5,5,0,1,0)
                g(1,3) = R_E(5,5,1,0,0)*R_E(5,5,0,0,1) + &
                    &Z_E(5,5,1,0,0)*Z_E(5,5,0,0,1)
                g(2,2) = R_E(5,5,0,1,0)*R_E(5,5,0,1,0) + &
                    &Z_E(5,5,0,1,0)*Z_E(5,5,0,1,0)
                g(2,3) = R_E(5,5,0,1,0)*R_E(5,5,0,0,1) + &
                    &Z_E(5,5,0,1,0)*Z_E(5,5,0,0,1)
                g(3,3) = R_E(5,5,0,0,1)*R_E(5,5,0,0,1) + &
                    &Z_E(5,5,0,0,1)*Z_E(5,5,0,0,1) + &
                    &R_E(5,5,0,0,0)*R_E(5,5,0,0,0)
                g(2,1) = g(1,2)
                g(3,1) = g(1,3)
                g(3,2) = g(2,3)
                call print_ar_2(g)
                write(*,*) 'rel diff'
                call print_ar_2(2*(g_E(5,5,:,:,0,0,0)-g)/(g_E(5,5,:,:,0,0,0)+g))
            end if
            
            do
                write(*,*) 'check D g_V?'
                if(yes_no(.false.)) then
                    do
                        write(*,*) 'derive in which coordinate? [1-3]'
                        read(*,*) id_d
                        if (id_d.ge.1 .and. id_d.le.3) then
                            exit
                        else
                            write(*,*) 'Choose between [1-3]'
                        end if
                    end do
                    d = 0
                    d(id_d) = 1
                    
                    write(*,*) 'D_'//trim(i2str(id_d))//' g_V(5,5) = '
                    call print_ar_2(g_E(5,5,:,:,d(1),d(2),d(3)))
                    
                    write(*,*) 'D_'//trim(i2str(id_d))//' g_V(5,5), manually'
                    g(1,1) = R_E(5,5,d(1)+1,d(2),d(3))*R_E(5,5,1,0,0) + &
                        &R_E(5,5,1,0,0)*R_E(5,5,d(1)+1,d(2),d(3)) + &
                        &Z_E(5,5,d(1)+1,d(2),d(3))*Z_E(5,5,1,0,0) + &
                        &Z_E(5,5,1,0,0)*Z_E(5,5,d(1)+1,d(2),d(3))
                    g(1,2) = R_E(5,5,d(1)+1,d(2),d(3))*R_E(5,5,0,1,0) + &
                        &R_E(5,5,1,0,0)*R_E(5,5,d(1),d(2)+1,d(3)) + &
                        &Z_E(5,5,d(1)+1,d(2),d(3))*Z_E(5,5,0,1,0) + &
                        &Z_E(5,5,1,0,0)*Z_E(5,5,d(1),d(2)+1,d(3))
                    g(1,3) = R_E(5,5,d(1)+1,d(2),d(3))*R_E(5,5,0,0,1) + &
                        &R_E(5,5,1,0,0)*R_E(5,5,d(1),d(2),d(3)+1) + &
                        &Z_E(5,5,d(1)+1,d(2),d(3))*Z_E(5,5,0,0,1) + &
                        &Z_E(5,5,1,0,0)*Z_E(5,5,d(1),d(2),d(3)+1)
                    g(2,2) = R_E(5,5,d(1),d(2)+1,d(3))*R_E(5,5,0,1,0) + &
                        &R_E(5,5,0,1,0)*R_E(5,5,d(1),d(2)+1,d(3)) + &
                        &Z_E(5,5,d(1),d(2)+1,d(3))*Z_E(5,5,0,1,0) + &
                        &Z_E(5,5,0,1,0)*Z_E(5,5,d(1),d(2)+1,d(3))
                    g(2,3) = R_E(5,5,d(1),d(2)+1,d(3))*R_E(5,5,0,0,1) + &
                        &R_E(5,5,0,1,0)*R_E(5,5,d(1),d(2),d(3)+1) + &
                        &Z_E(5,5,d(1),d(2)+1,d(3))*Z_E(5,5,0,0,1) + &
                        &Z_E(5,5,0,1,0)*Z_E(5,5,d(1),d(2),d(3)+1)
                    g(3,3) = R_E(5,5,d(1),d(2),d(3)+1)*R_E(5,5,0,0,1) + &
                        &R_E(5,5,0,0,1)*R_E(5,5,d(1),d(2),d(3)+1) + &
                        &Z_E(5,5,d(1),d(2),d(3)+1)*Z_E(5,5,0,0,1) + &
                        &Z_E(5,5,0,0,1)*Z_E(5,5,d(1),d(2),d(3)+1) + &
                        &2*R_E(5,5,d(1),d(2),d(3))*R_E(5,5,0,0,0)
                    g(2,1) = g(1,2)
                    g(3,1) = g(1,3)
                    g(3,2) = g(2,3)
                    call print_ar_2(g(:,:))
                    write(*,*) 'rel diff'
                    call print_ar_2(2*(g_E(5,5,:,:,d(1),d(2),d(3))-g)/&
                        &(g_E(5,5,:,:,d(1),d(2),d(3))+g))
                else
                    exit
                end if
            end do
                
            write(*,*) 'check g_F?'
            if(yes_no(.false.)) then
                write(*,*) 'g_F(5,5) = '
                call print_ar_2(g_F(5,5,:,:,0,0,0))
                
                write(*,*) 'THIS IS ONLY VALID FOR AXISYMMETRIC CASES'
                write(*,*) 'AND ONLY ELEMENT (1,3) IS VALID'
                write(*,*) 'g_F(5,5), manually'
                g = 0.0_dp
                !g(1,3) = -L_E(5,5,0,0,1)*(1-q_saf_E(5,0)*L_E(5,5,0,0,1))/&
                    !&(1+L_E(5,5,0,1,0))**2*g_E(5,5,2,2,0,0,0) + q_saf_E(5,0)* &
                    !&g_E(5,5,3,3,0,0,0) + (1-2*q_saf_E(5,0)*L_E(5,5,0,0,1))/&
                    !&(1+L_E(5,5,0,1,0))*g_E(5,5,2,3,0,0,0)
                g(1,3) = q_saf_E(5,0)*g_E(5,5,3,3,0,0,0)
                call print_ar_2(g)
                write(*,*) 'rel diff'
                call print_ar_2(2*(g_F(5,5,:,:,0,0,0)-g)/(g_F(5,5,:,:,0,0,0)+g))
            end if
            
            do
                write(*,*) 'check D g_F?'
                if(yes_no(.false.)) then
                    do
                        write(*,*) 'derive in which coordinate? [1-3]'
                        read(*,*) id_d
                        if (id_d.ge.1 .and. id_d.le.3) then
                            exit
                        else
                            write(*,*) 'Choose between [1-3]'
                        end if
                    end do
                    d = 0
                    d(id_d) = 1
                    
                    write(*,*) 'THIS IS ONLY VALID FOR AXISYMMETRIC CASES'
                    write(*,*) 'AND ONLY ELEMENT (1,3) IS VALID'
                    write(*,*) 'D_'//trim(i2str(id_d))//' g_F(5,5) = '
                    call print_ar_2(g_F(5,5,:,:,d(1),d(2),d(3)))
                    
                    write(*,*) 'D_'//trim(i2str(id_d))//' g_F(5,5), manually'
                    g = 0.0_dp
                    if (id_d.eq.1) then
                        g(1,3) = q_saf_E(5,1)*g_E(5,5,3,3,0,0,0) + &
                            &q_saf_E(5,0)*g_E(5,5,3,3,d(1),d(2),d(3))
                    else
                        g(1,3) = q_saf_E(5,0)*g_E(5,5,3,3,d(1),d(2),d(3))
                    end if
                    call print_ar_2(g)
                    write(*,*) 'rel diff'
                    call print_ar_2(2*(g_F(5,5,:,:,d(1),d(2),d(3))-g)/&
                        &(g_F(5,5,:,:,d(1),d(2),d(3))+g))
                else
                    exit
                end if
            end do
                
            
            write(*,*) 'Stopping'
            stop
        end if
    end function test_calc_g
    
    integer function test_calc_f_deriv() result(ierr)
        use metric_ops, only: calc_f_deriv
        use metric_vars, only: h_F, h_FD, T_FE
        use driver_rich, only: calc_eq
        use eq_vars, only: flux_p_E, flux_p_FD, grp_n_r_eq
        use X_vars, only: n_par
        
        character(*), parameter :: rout_name = 'test_calc_f_deriv'
        
        ! local variables
        real(dp) :: alt_calc(n_par, grp_n_r_eq)
        integer :: id, jd, kd
        integer :: der1(3)
        integer :: der2(3)
        integer :: der3(3)
        
        ! initialize ierr
        ierr = 0
        
        write(*,*) 'test calc_f_deriv?'
        if(yes_no(.false.)) then
            
            write(*,*) 'calculating equilibrium'
            ierr = calc_eq(0.2*pi)
            CHCKERR('')
            
            write(*,*) 'testing flux quantities'
            read(*,*)
            
            write(*,*) 'testing zeroth order'
            write(*,*) 'flux_p_FD = '
            call print_ar_1(flux_p_FD(:,0))
            write(*,*) 'alternative calculation = '
            alt_calc(1,:) = flux_p_E(:,0)
            call print_ar_1(alt_calc(1,:))
            write(*,*) 'rel difference = '
            call print_ar_1(2*abs(flux_p_FD(:,0)-alt_calc(1,:))/&
                &(flux_p_FD(:,0)+alt_calc(1,:)))
            read(*,*)
            
            write(*,*) 'testing first order'
            write(*,*) 'D_psi flux_p_FD = '
            call print_ar_1(flux_p_FD(:,1))
            write(*,*) 'alternative calculation = '
            alt_calc(1,:) = 2*pi
            call print_ar_1(alt_calc(1,:))
            write(*,*) 'rel difference = '
            call print_ar_1(2*abs(flux_p_FD(:,1)-alt_calc(1,:))/&
                &(flux_p_FD(:,1)+alt_calc(1,:)))
            read(*,*)
            
            write(*,*) 'testing second order'
            write(*,*) 'D^2_psi flux_p_FD = '
            call print_ar_1(flux_p_FD(:,2))
            write(*,*) 'alternative calculation = '
            alt_calc(1,:) = 0
            call print_ar_1(alt_calc(1,:))
            write(*,*) 'rel difference = '
            call print_ar_1(2*abs(flux_p_FD(:,2)-alt_calc(1,:))/&
                &(flux_p_FD(:,2)+alt_calc(1,:)))
            read(*,*)
            
            write(*,*) 'testing third order'
            write(*,*) 'D^3_psi flux_p_FD = '
            call print_ar_1(flux_p_FD(:,3))
            write(*,*) 'alternative calculation = '
            alt_calc(1,:) = 0
            call print_ar_1(alt_calc(1,:))
            write(*,*) 'rel difference = '
            call print_ar_1(2*abs(flux_p_FD(:,3)-alt_calc(1,:))/&
                &(flux_p_FD(:,3)+alt_calc(1,:)))
            read(*,*)
            
            
            write(*,*) 'testing normal quantities'
            read(*,*)
            
            write(*,*) 'testing zeroth order'
            write(*,*) 'h_FD(3,1) = '
            call print_ar_2(h_FD(:,:,3,1,0,0,0))
            write(*,*) 'alternative calculation = '
            alt_calc = h_F(:,:,3,1,0,0,0)
            call print_ar_2(alt_calc)
            write(*,*) 'rel difference = '
            call print_ar_2(2*abs(h_FD(:,:,3,1,0,0,0)-alt_calc)/&
                &(h_FD(:,:,3,1,0,0,0)+alt_calc))
            read(*,*)
            
            write(*,*) 'testing first order psi'
            write(*,*) 'D_psi h_FD(3,1) = '
            call print_ar_2(h_FD(:,:,3,1,0,1,0))
            write(*,*) 'alternative calculation = '
            alt_calc = 0.0_dp
            do id = 1,3
                der1 = 0
                der1(id) = 1
                alt_calc = alt_calc + T_FE(:,:,2,id,0,0,0)*&
                    &h_F(:,:,3,1,der1(1),der1(2),der1(3))
            end do
            call print_ar_2(alt_calc)
            write(*,*) 'rel difference = '
            call print_ar_2(2*abs(h_FD(:,:,3,1,0,1,0)-alt_calc)/&
                &(h_FD(:,:,3,1,0,1,0)+alt_calc))
            read(*,*)
            
            write(*,*) 'testing second order psi,psi'
            write(*,*) 'D_psi,psi h_FD(3,1) = '
            call print_ar_2(h_FD(:,:,3,1,0,2,0))
            write(*,*) 'alternative calculation = '
            alt_calc = 0.0_dp
            do id = 1,3
                do jd = 1,3
                    der1 = 0
                    der1(id) = 1
                    der2 = 0
                    der2(jd) = 1
                    alt_calc = alt_calc + T_FE(:,:,2,id,0,0,0)*&
                        &(T_FE(:,:,2,jd,der1(1),der1(2),der1(3))*&
                        &h_F(:,:,3,1,der2(1),der2(2),der2(3)) + &
                        &T_FE(:,:,2,jd,0,0,0)*h_F(:,:,3,1,der1(1)+der2(1),&
                        &der1(2)+der2(2),der1(3)+der2(3)))
                end do
            end do
            call print_ar_2(alt_calc)
            write(*,*) 'rel difference = '
            call print_ar_2(2*abs(h_FD(:,:,3,1,0,2,0)-alt_calc)/&
                &(h_FD(:,:,3,1,0,2,0)+alt_calc))
            read(*,*)
            
            write(*,*) 'testing second order psi,theta'
            write(*,*) 'D_psi,psi h_FD(3,1) = '
            call print_ar_2(h_FD(:,:,3,1,0,1,1))
            write(*,*) 'alternative calculation = '
            alt_calc = 0.0_dp
            do id = 1,3
                do jd = 1,3
                    der1 = 0
                    der1(id) = 1
                    der2 = 0
                    der2(jd) = 1
                    alt_calc = alt_calc + T_FE(:,:,3,id,0,0,0)*&
                        &(T_FE(:,:,2,jd,der1(1),der1(2),der1(3))*&
                        &h_F(:,:,3,1,der2(1),der2(2),der2(3)) + &
                        &T_FE(:,:,2,jd,0,0,0)*h_F(:,:,3,1,der1(1)+der2(1),&
                        &der1(2)+der2(2),der1(3)+der2(3)))
                end do
            end do
            call print_ar_2(alt_calc)
            write(*,*) 'rel difference = '
            call print_ar_2(2*abs(h_FD(:,:,3,1,0,1,1)-alt_calc)/&
                &(h_FD(:,:,3,1,0,1,1)+alt_calc))
            read(*,*)
            
            write(*,*) 'testing second order theta,psi'
            write(*,*) 'D_psi,psi h_FD(3,1) = '
            call print_ar_2(h_FD(:,:,3,1,0,1,1))
            write(*,*) 'alternative calculation = '
            alt_calc = 0.0_dp
            do id = 1,3
                do jd = 1,3
                    der1 = 0
                    der1(id) = 1
                    der2 = 0
                    der2(jd) = 1
                    alt_calc = alt_calc + T_FE(:,:,2,id,0,0,0)*&
                        &(T_FE(:,:,3,jd,der1(1),der1(2),der1(3))*&
                        &h_F(:,:,3,1,der2(1),der2(2),der2(3)) + &
                        &T_FE(:,:,3,jd,0,0,0)*h_F(:,:,3,1,der1(1)+der2(1),&
                        &der1(2)+der2(2),der1(3)+der2(3)))
                end do
            end do
            call print_ar_2(alt_calc)
            write(*,*) 'rel difference = '
            call print_ar_2(2*abs(h_FD(:,:,3,1,0,1,1)-alt_calc)/&
                &(h_FD(:,:,3,1,0,1,1)+alt_calc))
            read(*,*)
            
            write(*,*) 'testing third order psi,psi,theta'
            write(*,*) 'D_psi,psi,theta h_FD(3,1) = '
            call print_ar_2(h_FD(:,:,3,1,0,2,1))
            write(*,*) 'alternative calculation = '
            alt_calc = 0.0_dp
            do id = 1,3
                do jd = 1,3
                    do kd = 1,3
                        der1 = 0
                        der1(id) = 1
                        der2 = 0
                        der2(jd) = 1
                        der3 = 0
                        der3(kd) = 1
                        alt_calc = alt_calc + T_FE(:,:,3,id,0,0,0)*&
                            &(T_FE(:,:,2,jd,der1(1),der1(2),der1(3))*&
                            &  (T_FE(:,:,2,kd,der2(1),der2(2),der2(3))*&
                            &    h_F(:,:,3,1,der3(1),der3(2),der3(3)) + &
                            &    T_FE(:,:,2,kd,0,0,0)* &
                            &    h_F(:,:,3,1,der2(1)+der3(1),der2(2)+der3(2),der2(3)+der3(3)) &
                            &  ) + &
                            &  T_FE(:,:,2,jd,0,0,0) * &
                            &  (T_FE(:,:,2,kd,der1(1)+der2(1),der1(2)+der2(2),der1(3)+der2(3)) * &
                            &    h_F(:,:,3,1,der3(1),der3(2),der3(3)) + &
                            &    T_FE(:,:,2,kd,der2(1),der2(2),der2(3)) * &
                            &    h_F(:,:,3,1,der1(1)+der3(1),der1(2)+der3(2),der1(3)+der3(3)) + &
                            &    T_FE(:,:,2,kd,der1(1),der1(2),der1(3)) * &
                            &    h_F(:,:,3,1,der2(1)+der3(1),der2(2)+der3(2),der2(3)+der3(3)) + &
                            &    T_FE(:,:,2,kd,0,0,0) * &
                            &    h_F(:,:,3,1,der1(1)+der2(1)+der3(1),der1(2)+der2(2)+der3(2),der1(3)+der2(3)+der3(3)) &
                            &  )&
                            &)
                    end do
                end do
            end do
            call print_ar_2(alt_calc)
            write(*,*) 'rel difference = '
            call print_ar_2(2*abs(h_FD(:,:,3,1,0,2,1)-alt_calc)/&
                &(h_FD(:,:,3,1,0,2,1)+alt_calc))
            read(*,*)
            write(*,*) 'changing the order of the derivatives : '
            alt_calc = 0.0_dp
            do id = 1,3
                do jd = 1,3
                    do kd = 1,3
                        der1 = 0
                        der1(id) = 1
                        der2 = 0
                        der2(jd) = 1
                        der3 = 0
                        der3(kd) = 1
                        alt_calc = alt_calc + T_FE(:,:,2,id,0,0,0)*&
                            &(T_FE(:,:,3,jd,der1(1),der1(2),der1(3))*&
                            &  (T_FE(:,:,2,kd,der2(1),der2(2),der2(3))*&
                            &    h_F(:,:,3,1,der3(1),der3(2),der3(3)) + &
                            &    T_FE(:,:,2,kd,0,0,0)* &
                            &    h_F(:,:,3,1,der2(1)+der3(1),der2(2)+der3(2),der2(3)+der3(3)) &
                            &  ) + &
                            &  T_FE(:,:,3,jd,0,0,0) * &
                            &  (T_FE(:,:,2,kd,der1(1)+der2(1),der1(2)+der2(2),der1(3)+der2(3)) * &
                            &    h_F(:,:,3,1,der3(1),der3(2),der3(3)) + &
                            &    T_FE(:,:,2,kd,der2(1),der2(2),der2(3)) * &
                            &    h_F(:,:,3,1,der1(1)+der3(1),der1(2)+der3(2),der1(3)+der3(3)) + &
                            &    T_FE(:,:,2,kd,der1(1),der1(2),der1(3)) * &
                            &    h_F(:,:,3,1,der2(1)+der3(1),der2(2)+der3(2),der2(3)+der3(3)) + &
                            &    T_FE(:,:,2,kd,0,0,0) * &
                            &    h_F(:,:,3,1,der1(1)+der2(1)+der3(1),der1(2)+der2(2)+der3(2),der1(3)+der2(3)+der3(3)) &
                            &  )&
                            &)
                    end do
                end do
            end do
            call print_ar_2(alt_calc)
            write(*,*) 'rel difference = '
            call print_ar_2(2*abs(h_FD(:,:,3,1,0,2,1)-alt_calc)/&
                &(h_FD(:,:,3,1,0,2,1)+alt_calc))
            read(*,*)
            
            
            
            write(*,*) 'Stopping'
            stop
        end if
    end function test_calc_f_deriv
    
    integer function test_calc_inv_met() result(ierr)
        use metric_ops, only: calc_inv_met
        use metric_vars, only: T_EF, T_FE
        use driver_rich, only: calc_eq
        use eq_vars, only: grp_n_r_eq
        use X_vars, only: n_par
        
        character(*), parameter :: rout_name = 'test_calc_inv_met'
        
        ! local variables
        integer :: id, kd                                                       ! counters
        real(dp) :: TxT(3,3)
        real(dp) :: trace(n_par,grp_n_r_eq)
        
        ! initialize ierr
        ierr = 0
        
        write(*,*) 'test calc_inv_met?'
        if(yes_no(.false.)) then
            write(*,*) 'Testing whether T_VF is equal to T_FV'
            
            write(*,*) 'calculating equilibrium'
            ierr = calc_eq(0.2*pi)
            CHCKERR('')
            
            write(*,*) 'testing zeroth order'
            do kd = 1,grp_n_r_eq
                do id = 1,n_par
                    TxT = matmul(T_EF(id,kd,:,:,0,0,0),T_FE(id,kd,:,:,0,0,0))
                    trace(id,kd) = TxT(1,1)+TxT(2,2)+TxT(3,3)
                end do
            end do
            write(*,*) 'trace of T_VF*T_FV = '
            call print_ar_2(trace)
            read(*,*)
            
            write(*,*) 'testing first order: r derivative'
            do kd = 1,grp_n_r_eq
                do id = 1,n_par
                    TxT = matmul(T_EF(id,kd,:,:,1,0,0),T_FE(id,kd,:,:,0,0,0)) &
                        &+ matmul(T_EF(id,kd,:,:,0,0,0),T_FE(id,kd,:,:,1,0,0))
                    trace(id,kd) = sum(TxT)
                end do
            end do
            write(*,*) 'norm of Dr(T_VF*T_FV) = '
            call print_ar_2(trace)
            read(*,*)
            
            write(*,*) 'testing first order: t derivative'
            do kd = 1,grp_n_r_eq
                do id = 1,n_par
                    TxT = matmul(T_EF(id,kd,:,:,0,1,0),T_FE(id,kd,:,:,0,0,0)) &
                        &+ matmul(T_EF(id,kd,:,:,0,0,0),T_FE(id,kd,:,:,0,1,0))
                    trace(id,kd) = sum(TxT)
                end do
            end do
            write(*,*) 'norm of Dt(T_VF*T_FV) = '
            call print_ar_2(trace)
            read(*,*)
            
            write(*,*) 'testing first order: z derivative'
            do kd = 1,grp_n_r_eq
                do id = 1,n_par
                    TxT = matmul(T_EF(id,kd,:,:,0,0,1),T_FE(id,kd,:,:,0,0,0)) &
                        &+ matmul(T_EF(id,kd,:,:,0,0,0),T_FE(id,kd,:,:,0,0,1))
                    trace(id,kd) = sum(TxT)
                end do
            end do
            write(*,*) 'norm of Dz(T_VF*T_FV) = '
            call print_ar_2(trace)
            read(*,*)
            
            write(*,*) 'testing second order: tt derivative'
            do kd = 1,grp_n_r_eq
                do id = 1,n_par
                    TxT = matmul(T_EF(id,kd,:,:,0,2,0),T_FE(id,kd,:,:,0,0,0)) &
                        &+2*matmul(T_EF(id,kd,:,:,0,1,0),T_FE(id,kd,:,:,0,1,0))&
                        &+matmul(T_EF(id,kd,:,:,0,0,0),T_FE(id,kd,:,:,0,2,0))
                    trace(id,kd) = sum(TxT)
                end do
            end do
            write(*,*) 'norm of Dtt(T_VF*T_FV) = '
            call print_ar_2(trace)
            read(*,*)
            
            write(*,*) 'testing third order: rtt derivative'
            do kd = 1,grp_n_r_eq
                do id = 1,n_par
                    TxT = matmul(T_EF(id,kd,:,:,1,2,0),T_FE(id,kd,:,:,0,0,0)) &
                        &+matmul(T_EF(id,kd,:,:,0,2,0),T_FE(id,kd,:,:,1,0,0))&
                        &+2*matmul(T_EF(id,kd,:,:,1,1,0),T_FE(id,kd,:,:,0,1,0))&
                        &+2*matmul(T_EF(id,kd,:,:,0,1,0),T_FE(id,kd,:,:,1,1,0))&
                        &+matmul(T_EF(id,kd,:,:,1,0,0),T_FE(id,kd,:,:,0,2,0))&
                        &+matmul(T_EF(id,kd,:,:,0,0,0),T_FE(id,kd,:,:,1,2,0))
                    trace(id,kd) = sum(TxT)
                end do
            end do
            write(*,*) 'norm of Drtt(T_VF*T_FV) = '
            call print_ar_2(trace)
            read(*,*)
            
            
            write(*,*) 'Stopping'
            stop
        end if
    end function test_calc_inv_met

    integer function test_calc_det() result(ierr)
        use utilities, only: calc_det
        
        character(*), parameter :: rout_name = 'calc_det'
        
        ! local variables
        integer :: test_nr
        
        ! initialize ierr
        ierr = 0
        
        write(*,*) 'test calc_det?'
        if(yes_no(.false.)) then
            write(*,*) 'Testing whether calc_det is correctly calculated'
            do
                write(*,*) 'Select 1: 0D or 2: 2D'
                read(*,*) test_nr
                if (test_nr.gt.0 .and. test_nr.lt.3) then
                    exit
                else
                    write(*,*) 'You have to make a correct selection'
                end if
            end do
            
            select case (test_nr)
                case (1)
                    ierr = test_0()
                    CHCKERR('')
                case (2)
                    ierr = test_2()
                    CHCKERR('')
            end select
            
            
            write(*,*) 'Stopping'
            stop
        end if
    contains
        integer function test_0() result(ierr)
            character(*), parameter :: rout_name = 'test_0'
            
            ! local variables
            real(dp) :: A(3,3)
            real(dp) :: detA
            
            ! initialize ierr
            ierr = 0
            
            ! define a 3x3 matrix
            A = transpose(reshape([1,2,3,0,3,5,7,8,0],[3,3]))
            write(*,*) 'calculating determinant of:'
            call print_ar_2(A)
            
            detA = 0.0_dp
            ierr = calc_det(detA,A)
            CHCKERR('')
            write(*,*) 'det(A) = '//trim(r2str(detA))
            write(*,*) 'compare with analyitcal result: -33'
        end function test_0
        
        integer function test_2() result(ierr)
            character(*), parameter :: rout_name = 'test_2'
            
            ! local variables
            real(dp), allocatable :: A(:,:,:,:)
            real(dp), allocatable :: detA(:,:), detA_ALT(:,:)
            integer :: id, jd, kd                                               ! counters
            real(dp), allocatable :: A_part(:,:)
            
            ! initialize ierr
            ierr = 0
            
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
                    ierr = calc_det(detA(id,kd),A_part)
                    CHCKERR('')
                end do
            end do
            
            ! calculate alternative all at once 
            ierr = calc_det(detA_ALT,A)
            CHCKERR('')
            
            ! difference
            write(*,*) 'determinants calculated individually:'
            call print_ar_2(detA)
            write(*,*) 'determinants calculated global:'
            call print_ar_2(detA_ALT)
            write(*,*) 'difference'
            call print_ar_2(detA-detA_ALT)
        end function test_2
    end function test_calc_det
    
    integer function test_inv() result(ierr)
        use utilities, only: calc_inv
        
        character(*), parameter :: rout_name = 'test_inv'
        
        ! local variables
        real(dp), allocatable :: A(:,:,:,:), A_inv(:,:,:,:)
        integer :: id, jd, kd                                                   ! counters
        
        ! initialize ierr
        ierr = 0
        
        write(*,*) 'test inv?'
        if(yes_no(.false.)) then
            write(*,*) 'Testing whether inv works correctly'
            
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
            ierr = calc_inv(A_inv,A)
            CHCKERR('')
            
            ! multiply
            do kd = 1,5
                do id = 1,2
                    write(*,*) '(id,kd) = ('//trim(i2str(id))//','&
                        &//trim(i2str(kd))//')'
                    write(*,*) 'matrix A = '
                    call print_ar_2(A(id,kd,:,:))
                    write(*,*) 'inverse of A = '
                    call print_ar_2(A_inv(id,kd,:,:))
                    write(*,*) 'product = '
                    call print_ar_2(matmul(A(id,kd,:,:),A_inv(id,kd,:,:)))
                    write(*,*) ''
                    read(*,*)
                end do
            end do
            
            write(*,*) 'Stopping'
            stop
        end if
    end function test_inv
    
    integer function test_metric_transf() result(ierr)
        !use eq_ops, only: calc_eq
        !use X_vars, only: n_par
        !use eq_vars, only: flux_p_E, R_E, Z_E, L_E, theta, &
            !&zeta, &grp_n_r_eq
        !use VMEC, only: mpol, ntor, jac_V_c_H, jac_V_s_H, nfp, n_r_eq
        !use metric_ops, only: calc_inv_met
        !use metric_vars, only: jac_E, jac_F, T_FE, det_T_FE, g_F, h_F, g_E, T_EF
        !use utilities, only: calc_deriv, calc_det, conv_FHM
        !use num_vars, only: grid_style
        
        !character(*), parameter :: rout_name = 'test_metric_transf'
        
        !! local variables
        !real(dp) :: jac_ALT(1:n_par,1:grp_n_r_eq)                               ! VMEC calculated jac, FM
        !real(dp) :: jac_ALT_H(1:n_par,1:grp_n_r_eq)                             ! VMEC calculated jac, HM
        !real(dp) :: h_F_ALT(1:n_par,1:grp_n_r_eq)                               ! alternative of h_F (or derivatives)
        !integer :: id, jd, kd                                                   ! counter
        !real(dp) :: h_V(1:n_par,1:grp_n_r_eq,3,3,0:0,0:0,0:0)                   ! h_V
        !real(dp) :: gxh(1:n_par,1:grp_n_r_eq)                                   ! norm of g*h
        !real(dp) :: cs(0:mpol-1,-ntor:ntor,2)                                   ! (co)sines for all pol m and tor n
        !real(dp) :: zeta_in
        !real(dp), allocatable :: trigon_factors(:,:,:,:,:)                      ! trigonometric factor cosine for the inverse fourier transf.
        
        !! initialize ierr
        !ierr = 0
        
        !write(*,*) 'test metric_transf?'
        !if(yes_no(.false.)) then
            
            !output_i = 0
            
            !write(*,*) 'calculating equilibrium with constant zeta'
            !write(*,*) 'At which toroidal point zeta do you want the plot?'
            !read(*,*) zeta_in
            !grid_style = 2
            !ierr = calc_eq(zeta_in)
            !CHCKERR('')
            
            !write(*,*) 'check J_V?'
            !if(yes_no(.false.)) then
                !write(*,*) 'calculating jac_V directly to compare it with the 
                    !&calculation in the code')
                !! jac_V calculated
                !do kd = 1,grp_n_r_eq
                    !jac_ALT(:,kd) = R_E(:,kd,0,0,0)*(R_E(:,kd,1,0,0)*&
                        !&Z_E(:,kd,0,1,0)-R_E(:,kd,0,1,0)*&
                        !&Z_E(:,kd,1,0,0))
                !end do
                !call print_GP_3D('jac_V (1:calc, 2:direct calc, 3: diff)','',&
                    !&reshape([jac_E(:,:,0,0,0),jac_ALT,&
                    !&jac_E(:,:,0,0,0)-jac_ALT],[n_par,grp_n_r_eq,3]))
                
                !write(*,*) 'comparing jac_V with jac_V provided by VMEC (FM)'
                !! jac_V from VMEC
                !write(*,*) 'WRONG: YOU HAVE TO TAKE ONLY A PART OF THE NORMAL RANGE!!!'
                !ierr = calc_trigon_factors(theta,zeta,trigon_factors)
                !CHCKERR('')
                !ierr = fourier2real(jac_V_c_H,jac_V_s_H,trigon_factors,&
                    !&jac_ALT_H)
                !CHCKERR('')
                !! HM to FM
                !do id = 1,n_par
                    !ierr = conv_FHM(jac_ALT_H(id,:),jac_ALT(id,:),.false.)
                    !CHCKERR('')
                !end do
                !call print_GP_3D('jac_V (1: calc, 2: VMEC, 3: diff)','',&
                    !&reshape([jac_E(:,:,0,0,0),jac_ALT,&
                    !&jac_E(:,:,0,0,0)-jac_ALT],[n_par,grp_n_r_eq,3]))
                !call print_GP_3D('jac_V calc-VMEC [log]','',&
                    !&log10(max(2*abs(jac_E(:,:,0,0,0)-jac_ALT(:,:))/&
                    !&(jac_E(:,:,0,0,0)+jac_ALT),1E-10_dp)))
                
                !write(*,*) 'comparing jac_V with jac_V provided by VMEC (HM)'
                !! calculate half grid jac_V
                !do id = 1,n_par
                    !jac_ALT(id,2:grp_n_r_eq) = (n_r_eq-1._dp)*0.25*&
                        !&(R_E(id,1:grp_n_r_eq-1,0,0,0)+R_E(id,2:grp_n_r_eq,0,0,0))&
                        !&*((R_E(id,2:grp_n_r_eq,0,0,0)-R_E(id,1:grp_n_r_eq-1,0,0,0))&
                        !&  *(Z_E(id,1:grp_n_r_eq-1,0,1,0)+Z_E(id,2:grp_n_r_eq,0,1,0)) &
                        !&-(Z_E(id,2:grp_n_r_eq,0,0,0)-Z_E(id,1:grp_n_r_eq-1,0,0,0))&
                        !&  *(R_E(id,1:grp_n_r_eq-1,0,1,0)+R_E(id,2:grp_n_r_eq,0,1,0)))
                !end do
                !jac_ALT(:,1) = 0.0_dp
                !call print_GP_3D('(HM) jac_V (1: calc, 2: VMEC, 3: diff)','',&
                    !&reshape([jac_ALT,jac_ALT_H,&
                    !&jac_ALT-jac_ALT_H],[n_par,grp_n_r_eq,3]))
                !call print_GP_3D('(HM) jac_V calc-VMEC [log]','',&
                    !&log10(max(2*abs(jac_ALT-jac_ALT_H(:,:))/&
                    !&(jac_ALT+jac_ALT_H),1E-10_dp)))
                !call print_GP_3D('(HM) jac_V calc-VMEC','',&
                    !&jac_ALT-jac_ALT_H)
                !call print_GP_2D('(HM) jac_V calc-VMEC','',&
                    !&transpose(jac_ALT-jac_ALT_H))
            !end if
            
            !write(*,*) 'check J_F?'
            !if(yes_no(.false.)) then
                !!   jac_F = (flux_p_E'/2pi * (1+L_t))^-1 * R * (R'Z_t - R' Z_t)
                !do kd = 1,grp_n_r_eq
                    !jac_ALT(:,kd) = R_E(:,kd,0,0,0)*(R_E(:,kd,1,0,0)*&
                        !&Z_E(:,kd,0,1,0)-R_E(:,kd,0,1,0)*&
                        !&Z_E(:,kd,1,0,0))/(flux_p_E(kd,1)/(2*pi)*&
                        !&(1+L_E(:,kd,0,1,0)))
                !end do
                !call print_GP_3D('jac_F (1:calc, 2:direct calc, 3: diff)','',&
                    !&reshape([jac_F(:,:,0,0,0),jac_ALT,&
                    !&jac_F(:,:,0,0,0)-jac_ALT],[n_par,grp_n_r_eq,3]))
            !end if
            
            !!write(*,*) 'check r and rr derivatives?'
            !!if(yes_no(.false.)) then
                !!call print_GP_2D('J_F(par-'//trim(i2str(par))//')','',&
                    !!&jac_F(par,:,0,0,0))
                !!do kd = 1,2
                    !!write(*,*) 'checking r^'//trim(i2str(kd))//' deriv.'
                    !!ierr = calc_deriv(jac_F(par,:,0,0,0),&
                        !!&jac_ALT(par,:),n_r-1._dp,kd,1)
                    !!CHCKERR('')
                    !!call print_GP_2D('Dr^'//trim(i2str(kd))//' J_F(par=par) &
                        !!&(1: an, 2: num)','',&
                        !!&reshape([jac_F(par,:,kd,0,0),jac_ALT(par,:)],[grp_n_r_eq,2]))
                    !!call print_GP_2D('diff with num Dr^'//trim(i2str(kd))//' &
                        !!&J_F(par='//trim(i2str(par))//') [log]','',&
                        !!&log10(max(2*abs(jac_ALT(par,:)-&
                        !!&jac_F(par,:,kd,0,0))/(jac_ALT(par,:)+&
                        !!&jac_F(par,:,kd,0,0)),1E-10_dp)))
                !!end do
            !!end if
            
            !write(*,*) 'check jacobians as determinants?'
            !if(yes_no(.false.)) then
                !ierr = calc_det(jac_ALT,T_FE(:,:,:,:,0,0,0))
                !CHCKERR('')
                !call print_ar_2(jac_ALT-det_T_FE(:,:,0,0,0))
                !read(*,*)
            !end if
            
            !write(*,*) 'check g*h?'
            !if(yes_no(.false.)) then
                !write(*,*) 'VMEC coordinate system'
                !! calculate h_V from g_E
                !ierr = calc_inv_met(h_V,g_E,[0,0,0])
                !CHCKERR('')
                !do kd = 1,grp_n_r_eq
                    !do id = 1,n_par
                        !!write(*,*) 'h_V('//trim(i2str(id))//','//
                            !!&trim(i2str(kd))//') = ')
                        !!call print_ar_2(h_V(id,kd,:,:,0,0,0))
                        !!write(*,*) 'g_V('//trim(i2str(id))//','//
                            !!&trim(i2str(kd))//') = ')
                        !!call print_ar_2(g_E(id,kd,:,:,0,0,0))
                        !!write(*,*) 'h_V*g_V'
                        !!call print_ar_2(matmul(h_V(id,kd,:,:,0,0,0),&
                            !!&g_E(id,kd,:,:,0,0,0)))
                        !gxh(id,kd) = sum(matmul(h_V(id,kd,:,:,0,0,0),&
                            !&g_E(id,kd,:,:,0,0,0))-&
                            !&reshape([1,0,0,0,1,0,0,0,1],[3,3]))
                    !end do
                !end do
                !write(*,*) 'g_V * h_V = '
                !call print_ar_2(gxh)
                !read(*,*)
                
                !write(*,*) 'flux coordinate system'
                !do kd = 1,grp_n_r_eq
                    !do id = 1,n_par
                        !gxh(id,kd) = sum(matmul(h_F(id,kd,:,:,0,0,0),&
                            !&g_F(id,kd,:,:,0,0,0))-&
                            !&reshape([1,0,0,0,1,0,0,0,1],[3,3]))
                    !end do
                !end do
                !write(*,*) 'g_F * h_F = '
                !call print_ar_2(gxh)
                !read(*,*)
            !end if
            
            !write(*,*) 'check T_XY*T_YX?'
            !if(yes_no(.false.)) then
                !write(*,*) 'flux coordinate system'
                !do kd = 1,grp_n_r_eq
                    !do id = 1,n_par
                        !gxh(id,kd) = sum(matmul(T_FE(id,kd,:,:,0,0,0),&
                            !&T_EF(id,kd,:,:,0,0,0))-&
                            !&reshape([1,0,0,0,1,0,0,0,1],[3,3]))
                    !end do
                !end do
                !write(*,*) 'T_FV * T_VF = '
                !call print_ar_2(gxh)
                !read(*,*)
            !end if
            
            !write(*,*) 'check h_F?'
            !if(yes_no(.false.)) then
                !do kd = 1,3
                    !do id = 1,3
                        !write(*,*) 'checking r derivatives'
                        !call print_GP_3D('h_F('//trim(i2str(id))//','//&
                            !&trim(i2str(kd))//')','',h_F(:,2:grp_n_r_eq,id,kd,0,0,0))
                        !do jd =1,n_par
                            !ierr = calc_deriv(h_F(jd,2:grp_n_r_eq,id,kd,0,0,0),&
                                !&h_F_ALT(jd,2:grp_n_r_eq),n_r_eq-1._dp,1,1)
                            !CHCKERR('')
                        !end do
                        !call print_GP_3D('Dr h_F('//trim(i2str(id))//','//&
                            !&trim(i2str(kd))//') (1: an, 2: num, 3:diff)',&
                            !&'',reshape([h_F(:,2:grp_n_r_eq,id,kd,1,0,0),&
                            !&h_F_ALT(:,2:grp_n_r_eq),h_F(:,2:grp_n_r_eq,id,kd,1,0,0)-&
                            !&h_F_ALT(:,2:grp_n_r_eq)],[n_par,grp_n_r_eq-1,3]))
                        !call print_GP_3D('Dr h_F an-num [log]','',log10(&
                            !&max(abs(2*(h_F_ALT(:,2:grp_n_r_eq)-&
                            !&h_F(:,2:grp_n_r_eq,id,kd,1,0,0))/(maxval(h_F_ALT(:,2:grp_n_r_eq)+&
                            !&h_F(:,2:grp_n_r_eq,id,kd,1,0,0)))),1E-10_dp)))
                    !end do
                !end do
            !end if
                
            
            !write(*,*) 'Stopping'
            !stop
        !end if
    end function test_metric_transf

    integer function test_B() result(ierr)
        use X_vars, only: n_par
        use eq_vars, only: q_saf_E, flux_p_E, L_E, theta_E, zeta_E, &
            &grp_n_r_eq, grp_min_r_eq, grp_max_r_eq
        !use eq_vars, only: theta, zeta, pres_FD
        use VMEC, only: B_V_sub_c_M, B_V_sub_s_M, B_V_s_H, B_V_c_H
        !use VMEC, only: mpol, ntor, B_V_sub_c_M, &
            !&B_V_sub_s_M, nfp
        use metric_vars, only: g_E, jac_E, g_F, jac_F, T_FE
        use metric_ops, only: calc_inv_met
        !use metric_vars, only: g_FD, jac_FD
        use utilities, only: calc_deriv, conv_FHM
        use driver_rich, only: calc_eq
        use num_vars, only: grid_style
        use MPI_ops, only: split_MPI
        use fourier_ops, only: calc_trigon_factors, fourier2real
        
        character(*), parameter :: rout_name = 'test_B'
        
        ! local variables
        integer :: id, jd, kd                                                   ! counters
        real(dp), allocatable :: B(:,:)                                         ! magn. field
        real(dp), allocatable :: B_ALT(:,:)                                     ! magn. field alternative
        real(dp), allocatable :: B_ALT_2(:,:)                                   ! magn. field alternative 2
        real(dp), allocatable :: B_V_sub(:,:,:)                                 ! cov. components of B in VMEC coords
        real(dp), allocatable :: B_V_sub_ALT(:,:,:)                             ! VMEC output
        real(dp), allocatable :: B_F_sub(:,:,:)                                 ! cov. components of B in flux coords
        real(dp), allocatable :: B_F_sub_ALT(:,:,:)                             ! conversion of VMEC output
        real(dp), allocatable :: h_V(:,:,:,:,:,:,:)                             ! h_V
        real(dp), allocatable :: dum(:,:)
        real(dp), allocatable :: dum2(:,:)
        real(dp), allocatable :: dum3(:,:)
        real(dp), allocatable :: trigon_factors(:,:,:,:,:)                      ! trigonometric factor cosine for the inverse fourier transf.
        
        ! initialize ierr
        ierr = 0
        
        write(*,*) 'test B?'
        if(yes_no(.false.)) then
            
            write(*,*) 'calculating equilibrium with constant zeta'
            grid_style = 2
            ierr = split_MPI()
            CHCKERR('')
            ierr = calc_eq(0.12*pi)
            CHCKERR('')
            
            ! allocate variables
            allocate(B(1:n_par,1:grp_n_r_eq))
            allocate(B_ALT(1:n_par,1:grp_n_r_eq))
            allocate(B_ALT_2(1:n_par,1:grp_n_r_eq))
            allocate(B_V_sub(1:n_par,1:grp_n_r_eq,3))
            allocate(B_V_sub_ALT(1:n_par,1:grp_n_r_eq,3))
            allocate(B_F_sub(1:n_par,1:grp_n_r_eq,3))
            allocate(B_F_sub_ALT(1:n_par,1:grp_n_r_eq,3))
            allocate(h_V(1:n_par,1:grp_n_r_eq,3,3,0:0,0:0,0:0))
            allocate(dum(1:n_par,1:grp_n_r_eq))
            allocate(dum2(1:n_par,1:grp_n_r_eq))
            allocate(dum3(1:n_par,1:grp_n_r_eq))
            
            write(*,*) 'check covar. comp. of magnetic field in VMEC coords?'
            if(yes_no(.false.)) then
                ierr = calc_B_V_covar([0,0])
                CHCKERR('')
                do jd = 1,3
                    call print_GP_2D('B_'//trim(i2str(jd))//'(par=5) &
                        &(1: calc, 2: VMEC)','',&
                        &reshape([B_V_sub(5,:,jd),B_V_sub_ALT(5,:,jd)],[grp_n_r_eq,2]))
                    call print_GP_2D('B_'//trim(i2str(jd))//' - B_ALT_'//&
                        &trim(i2str(jd))//'(par=5) from VMEC, REL error [log]',&
                        &'',log10(max(abs(2*(B_V_sub(5,2:grp_n_r_eq,jd)-&
                        &B_V_sub_ALT(5,2:grp_n_r_eq,jd))/(B_V_sub(5,2:grp_n_r_eq,jd)+&
                        &B_V_sub_ALT(5,2:grp_n_r_eq,jd))),1E-10_dp)))
                end do
            end if
            
            write(*,*) 'check covar. comp. of magnetic field in flux coords?'
            if(yes_no(.false.)) then
                write(*,*) 'calculating B_F_sub from g_theta,i'
                do jd = 1,3
                    B_F_sub(:,:,jd) = g_F(:,:,3,jd,0,0,0)/jac_F(:,:,0,0,0)
                end do
                
                write(*,*) 'calculating B_F_sub_ALT from B_V_sub_ALT from VMEC'
                ierr = calc_B_V_covar([0,0])                                    ! calculating B_V_sub_ALT
                CHCKERR('')
                do jd = 1,3
                    B_F_sub_ALT(:,:,jd) = 0.0_dp                                !  transforming B_V_sub_ALT into B_F_sub_ALT
                    do id = 1,3
                        B_F_sub_ALT(:,:,jd) = B_F_sub_ALT(:,:,jd) + &
                            &T_FE(:,:,jd,id,0,0,0)*B_V_sub_ALT(:,:,id)
                    end do
                end do
                
                do jd = 1,3
                    call print_GP_2D('B_'//trim(i2str(jd))//'(par=5) &
                        &(1: calc from g_F, 2: VMEC)','',&
                        &reshape([B_F_sub(5,:,jd),B_F_sub_ALT(5,:,jd)],[grp_n_r_eq,2]))
                    call print_GP_2D('B_'//trim(i2str(jd))//' - B_ALT_'//&
                        &trim(i2str(jd))//'(par=5) from VMEC, REL error [log]',&
                        &'',log10(max(abs(2*(B_F_sub(5,2:grp_n_r_eq,jd)-&
                        &B_F_sub_ALT(5,2:grp_n_r_eq,jd))/(B_F_sub(5,2:grp_n_r_eq,jd)+&
                        &B_F_sub_ALT(5,2:grp_n_r_eq,jd))),1E-10_dp)))
                end do
            end if
            
            write(*,*) 'check magn. of. magnetic field?'
            if(yes_no(.false.)) then
                write(*,*) 'comparison with g_V(1,2)/J^2'
                ierr = calc_trigon_factors(theta_E,zeta_E,trigon_factors)
                CHCKERR('')
                ierr = fourier2real(&
                    &B_V_c_H(:,:,grp_min_r_eq:grp_max_r_eq), &
                    &B_V_s_H(:,:,grp_min_r_eq:grp_max_r_eq), &
                    &trigon_factors,B_ALT(:,:))
                CHCKERR('')
                deallocate(trigon_factors)
                B(:,:) = sqrt(g_F(:,:,3,3,0,0,0))/&
                    &jac_F(:,:,0,0,0)
                call print_GP_2D('B(par=5) &
                    &(1: calc from g_F, 2: VMEC)','',reshape([B(5,:),&
                    &B_ALT(5,:)],[grp_n_r_eq,2]))
                call print_GP_2D('B - B_ALT (par=5) from VMEC, REL error [log]'&
                    &,'',log10(max(abs(2*(B(5,2:grp_n_r_eq)-B_ALT(5,2:grp_n_r_eq))/&
                    &(B(5,2:grp_n_r_eq)+B_ALT(5,2:grp_n_r_eq))),1E-10_dp)))
                
                write(*,*) 'comparison with B_i B_j h^ij'
                ! calculate h_V from g_V
                ierr = calc_inv_met(h_V,g_E,[0,0,0])
                CHCKERR('')
                ierr = calc_B_V_covar([0,0])
                CHCKERR('')
                B_ALT_2 = 0.0_dp
                do id = 1,3
                    do kd = 1,3
                        B_ALT_2(:,:) = B_ALT_2(:,:) + &
                            &B_V_sub(:,:,id)*B_V_sub(:,:,kd)*&
                            &h_V(:,:,id,kd,0,0,0)
                    end do
                end do
                B_ALT_2 = sqrt(B_ALT_2)
                ! conversion FM -> HM and store in B
                do id = 1,n_par
                    ierr = conv_FHM(B_ALT_2(id,:),B(id,:),.true.)
                    CHCKERR('')
                end do
                call print_GP_2D('B(par=5) (1: calc from B_i B^j, 2: VMEC)',&
                    &'',reshape([B(5,:),B_ALT(5,:)],[grp_n_r_eq,2]))
                call print_GP_2D('B - B_ALT (par=5) from VMEC, REL error [log]'&
                    &,'',log10(max(abs(2*(B(5,2:grp_n_r_eq)-B_ALT(5,2:grp_n_r_eq))/&
                    &(B(5,2:grp_n_r_eq)+B_ALT(5,2:grp_n_r_eq))),1E-10_dp)))
            end if
            
            write(*,*) 'test whether D_alpha B_theta = D_theta B_alpha?'
            if(yes_no(.false.)) then
                
                write(*,*) 'plot B_zeta'
                dum(:,5) = flux_p_E(5,1)/(2*pi*jac_E(:,5,0,0,0)) * &
                    &(-q_saf_E(5,0)*(1+L_E(:,5,0,1,0))*g_E(:,5,3,3,0,0,0) - &
                    &(1-q_saf_E(5,0)*L_E(:,5,0,0,1))*g_E(:,5,3,2,0,0,0))
                write(*,*) 'theta (should be increasing in the first dim.) = '
                call print_ar_2(theta_E)
                write(*,*) 'zeta (should be constant) = '
                call print_ar_2(zeta_E)
                write(*,*) 'B_zeta from VMEC = '
                ierr = calc_B_V_covar([0,0])                                    ! calculating B_V_sub_ALT
                CHCKERR('')
                call print_ar_2(B_V_sub_ALT(:,:,3))
                call print_GP_2D('B_zeta (1: calc, 2: VMEC)','',&
                    &reshape([dum(:,5),B_V_sub_ALT(:,5,3)],[n_par,2]))
                call print_GP_2D('B_zeta difference, calc-VMEC','',&
                    &dum(:,5)-B_V_sub_ALT(:,5,3))
            
                !write(*,*) 'plot D_theta B_zeta'
                !dum(:,5) = flux_p_E(5,1)/(2*pi*jac_E(:,5,0,0,0)) * &
                    !&(-jac_E(:,5,0,1,0)/jac_E(:,5,0,0,0) * &
                    !&(q_saf_E(5,0)*(1+L_E(:,5,0,1,0))*g_E(:,5,3,3,0,0,0) + &
                    !&(1-q_saf_E(5,0)*L_E(:,5,0,0,1))*g_E(:,5,3,2,0,0,0)) + &
                    !&q_saf_E(5,0)*(L_E(:,5,0,2,0)*g_E(:,5,3,3,0,0,0)+&
                    !&(1+L_E(:,5,0,1,0))*g_E(:,5,3,3,0,1,0)) &
                    !&-q_saf_E(5,0)*L_E(:,5,0,1,1)*g_E(:,5,3,2,0,0,0) + &
                    !&(1-q_saf_E(5,0)*L_E(:,5,0,0,1))*g_E(:,5,3,2,0,1,0))
                !ierr = calc_B_V_covar([1,0])
                !CHCKERR('')
                !call print_GP_2D('D_theta B_zeta (1: calc, 2: an)','',&
                    !&reshape([dum(:,5),B_V_sub_ALT(:,5,3)],[n_par,2]))
                !call print_GP_2D('D_theta B_zeta difference, calc-VMEC','',&
                    !&dum(:,5)-B_V_sub_ALT(:,5,3))
                !call print_GP_2D('D_theta B_zeta REL difference, calc-VMEC','',&
                    !&abs(dum(:,5)-B_V_sub_ALT(:,5,3))/maxval(dum(:,5)))
                
                !write(*,*) 'D_alpha B_theta = '
                !dum = g_FD(:,:,3,3,1,0,0)/jac_FD(:,:,0,0,0)- &
                    !&g_FD(:,:,3,3,0,0,0)*jac_FD(:,:,1,0,0)/&
                    !&jac_FD(:,:,0,0,0)**2
                !call print_ar_2(dum)
                !write(*,*) 'D_theta B_alpha = '
                !dum2 = g_FD(:,:,3,1,0,0,1)/jac_FD(:,:,0,0,0)- &
                    !&g_FD(:,:,3,1,0,0,0)*jac_FD(:,:,0,0,1)/&
                    !&jac_FD(:,:,0,0,0)**2
                !call print_ar_2(dum2)
                !write(*,*) 'D_alpha B_theta - D_theta B_alpha 
                    !&[rel. error] = ')
                !call print_ar_2(2*abs((dum-dum2)/(dum+dum2)))
                !call print_GP_2D('rel error at par = 5, [log]','',&
                    !&log10(max(abs((dum(5,:)-dum2(5,:))/(dum(5,:)+&
                    !&dum2(5,:))),1E-10_dp)))
                !!write(*,*) 'D_theta B_alpha converted to VMEC coords = '
                !!ierr = calc_B_V_covar([1,0])                                    ! so B_V_sub_ALT contains the theta derivatives
                !!CHCKERR('')
                !!call print_ar_2(1/(1+L_E(:,:,0,1,0))*B_V_sub_ALT(:,:,3))
                !read(*,*)
                
            end if
            
            !write(*,*) 'test whether D_theta B_psi - D_psi B_theta = mu_0 J &
                !&p''?'
            !if(yes_no(.false.)) then
            !write(*,*) 'ADAPT THIS!!!!! MU_0 GOES AWAY FOR NORMALIZED QUANTITIES !!!!!!!!!!'
                    !! dum = D_theta B_psi
                    !dum = g_FD(:,:,3,2,0,0,1)/jac_FD(:,:,0,0,0)- &
                        !&g_FD(:,:,3,2,0,0,0)*jac_FD(:,:,0,0,1)/&
                        !&jac_FD(:,:,0,0,0)**2
                    !dum(:,1) = 0.0_dp
                    !! dum2 = D_psi B_theta
                    !dum2 = g_FD(:,:,3,3,0,1,0)/jac_FD(:,:,0,0,0)- &
                        !&g_FD(:,:,3,3,0,0,0)*jac_FD(:,:,0,1,0)/&
                        !&jac_FD(:,:,0,0,0)**2
                    !dum2(:,1) = 0.0_dp
                    !call print_GP_2D('1: D_theta B_psi, 2: D_psi B_theta','',&
                        !&reshape([dum(5,:),dum2(5,:)],[grp_n_r_eq,2]))
                    
                    !do id = 1,n_par
                        !dum3(id,:) = mu_0*jac_FD(id,:,0,0,0)*pres_FD(:,1)
                    !end do
                    !call print_GP_2D('1: D_theta B_psi-D_psi B_theta, 2: calc',&
                        !&'',reshape([dum(5,:)-dum2(5,:),dum3(5,:)],[grp_n_r_eq,2]))
                    !call print_GP_2D('diff','',dum(5,:)-dum2(5,:)-dum3(5,:))
            !end if
            
            
            write(*,*) 'Stopping'
            stop
        end if
    contains
        ! calculates B_V_sub  from VMEC output and also  B_V_sub_ALT from direct
        ! calculation using the metric factors
        integer function calc_B_V_covar(deriv) result(ierr)
            character(*), parameter :: rout_name = 'calc_B_V_covar'
            
            ! input / output
            integer :: deriv(2)
            
            ! local variables
            real(dp) :: tempvar(1:n_par,1:grp_n_r_eq)                           ! for HM to FM conversions
            
            ! initialize ierr
            ierr = 0
            
            if (deriv(1).ne.0 .or. deriv(2).ne.0) write(*,*) 'WARNING: In calc&
                &_B_V_covar, the returned value of "B_V_sub" is NOT correct &
                &because derivatives are asked'
            
            ierr = calc_trigon_factors(theta_E,zeta_E,trigon_factors)
            CHCKERR('')
            do jd = 1,3
                ierr = fourier2real(&
                    &B_V_sub_c_M(:,:,grp_min_r_eq:grp_max_r_eq,jd), &
                    &B_V_sub_s_M(:,:,grp_min_r_eq:grp_max_r_eq,jd), &
                    &trigon_factors,B_V_sub_ALT(:,:,jd),deriv)
                CHCKERR('')
                do kd = 1,grp_n_r_eq
                    B_V_sub(:,kd,jd) = flux_p_E(kd,1)/&
                        &(2*pi*jac_E(:,kd,0,0,0)) * &
                        &(-q_saf_E(kd,0)*(1+L_E(:,kd,0,1,0)) * &
                        &g_E(:,kd,3,jd,0,0,0) - (1-q_saf_E(kd,0)*&
                        &L_E(:,kd,0,0,1))*g_E(:,kd,2,jd,0,0,0))
                end do
            end do
            deallocate(trigon_factors)
            
            ! convert HM to FM (only theta and zeta component)
            do jd = 2,3
                do id = 1,n_par
                    ierr = conv_FHM(B_V_sub_ALT(id,:,jd),tempvar(id,:),&
                        &.false.)
                    CHCKERR('')
                end do
                B_V_sub_ALT(:,:,jd) = tempvar
            end do
        end function calc_B_V_covar
    end function test_B
    
    integer function test_pres_balance() result(ierr)
        use MPI_ops, only: split_MPI
        use driver_rich, only: calc_eq
        use X_vars, only: n_par
        use eq_vars, only: grp_n_r_eq
        
        character(*), parameter :: rout_name = 'test_pres_balance'
        
        ! local variables
        real(dp), allocatable :: B(:,:,:,:,:,:)                                 ! magnetic field
        
        ! initialize ierr
        ierr = 0
        
        write(*,*) 'test pressure balance?'
        if(yes_no(.false.)) then
            
            write(*,*) 'testing pressure balance'
            
            write(*,*) 'calculating metric coefficients'
            
            ierr = split_MPI()
            CHCKERR('')
            ierr = calc_eq(0.12*pi)
            CHCKERR('')
            
            write(*,*) 'calculating magnetic field and derivatives'
            allocate(B(n_par,grp_n_r_eq,0:1,0:1,0:1,3))
        end if
    end function test_pres_balance

    integer function test_ang_B() result(ierr)
        !use VMEC, only: mpol, ntor
        !use X_vars, only: n_par
        !use eq_vars, only: calc_eqd_grid, calc_ang_grid, &
            !&theta, zeta, grp_n_r_eq
        !use num_vars, only: use_tor_flux
        
        !character(*), parameter :: rout_name = 'test_ang_B'
        
        !! local variables (not to be used in child functions)
        !real(dp) :: alpha
        !integer :: id
        !real(dp), allocatable :: plot_dep_var(:)
        !real(dp), allocatable :: plot_var(:,:)
        !integer :: n_plot
        !character(len=max_str_ln) :: par_ang, dep_ang                           ! parallel angle and dependent angle, for output message
        !real(dp) :: grid_min, grid_max
        
        !! local variables (also used in child functions)
        !integer :: kd
        
        !! initialize ierr
        !ierr = 0
        
        !write(*,*) 'test ang_B?'
        !if(yes_no(.false.)) then
            !output_i = 0
            !write(*,*) 'Plot zeta(theta)'
            
            
            !! CALCULATE THETA(ZETA)
            !write(*,*) 'Calculating theta(zeta)'
            !alpha = 1.2_dp*pi
            
            !! decrease the number of parallel points
            !n_par = 100
            
            !! calculate grid points (theta, zeta) that follow the magnetic field
            !! line
            !ierr = calc_ang_grid(alpha)
            !CHCKERR('')
            
            !call print_GP_2D('zeta(theta)','',zeta,x=theta)
            !!do kd = 1,grp_n_r_eq
                !!call print_GP_2D('zeta(theta) at r = '//trim(i2str(kd))//'/'//&
                    !!&trim(i2str(grp_n_r_eq)),'',zeta(:,kd),x=theta(:,kd))
            !!end do
            
            
            !! CALCULATE F FOR A RANGE OF PARALLEL VALUES
            !write(*,*) 'Calculating f for a range of parallel values'
            
            !n_plot = 100
            
            !! initialize
            !allocate(plot_dep_var(n_plot)); plot_dep_var = 0.0_dp
            !allocate(plot_var(2,n_plot)); plot_var = 0.0_dp
            
            !! set correct plot messages
            !if (use_tor_flux) then                                              ! looking for zeta
                !par_ang = 'theta'; dep_ang = 'zeta'
            !else                                                                ! looking for theta
                !par_ang = 'zeta'; dep_ang = 'theta'
            !end if
                
            !do
                !write(*,*) 'At which grp_n_r_eq do you want to plot the solutions? [1..'&
                    !&//trim(i2str(grp_n_r_eq))//']'
                !read(*,*) kd
                !if (kd.ge.1 .and. kd.le.grp_n_r_eq) then
                    !exit
                !else
                    !write(*,*) 'Please choose an acceptable value'
                !end if
            !end do

            !! determine grid considering the solutions on current flux surface
            !if (use_tor_flux) then                                              ! looking for zeta
                !grid_min = minval(zeta(:,kd)) - &
                    !&0.1*(maxval(zeta(:,kd))-minval(zeta(:,kd)))
                !grid_max = maxval(zeta(:,kd)) + &
                    !&0.1*(maxval(zeta(:,kd))-minval(zeta(:,kd)))
            !else                                                                ! looking for theta
                !grid_min = minval(theta(:,kd)) - &
                    !&0.1*(maxval(theta(:,kd))-minval(theta(:,kd)))
                !grid_max = maxval(theta(:,kd)) + &
                    !&0.1*(maxval(theta(:,kd))-minval(theta(:,kd)))
            !end if
            !ierr = calc_eqd_grid(plot_dep_var,n_plot,grid_min,grid_max)
            !CHCKERR('')
            !par: do id = 1, n_par
                !if (use_tor_flux) then                                          ! looking for zeta
                    !plot_var(1,:) = plot_dep_var                                ! eq. grid of zeta
                    !plot_var(2,:) = find_f_plot(n_plot,plot_dep_var, &
                        !&theta(id,kd),alpha)                                    ! corresponding f(zeta)
                !else                                                            ! looking for theta
                    !plot_var(1,:) = plot_dep_var                                ! eq. grid of theta
                    !plot_var(2,:) = find_f_plot(n_plot,plot_dep_var, &
                        !&zeta(id,kd),alpha)                                     ! corresponding f(theta)
                !end if
                
                !call print_GP_2D('zeta - q (theta+lambda) - alpha_0 at &
                    !&(r,theta,zeta) = ('//&
                    !&trim(r2strt((kd-1._dp)/(grp_n_r_eq-1._dp)))//','//&
                    !&trim(r2strt(theta(id,kd)))//', '//&
                    !&trim(r2strt(zeta(id,kd)))//')','',plot_var(2,:),&
                    !&x=plot_var(1,:))
                
                !write(*,*) 'Paused... plot next?'
                !if(.not.yes_no(.true.)) then
                    !exit
                !end if
            !end do par
            
            
            !write(*,*) 'Stopping'
            !stop
        !end if
    !contains
        !! calculates the function f = zeta - q (theta + lambda) - alpha_0
        !! makes use of kd from parent to indicate the flux surface
        !function find_f_plot(n_ang,dep_var,par_var,alpha)
            !use fourier_ops, only: calc_ang_grid, f2r
            !use VMEC, only: iotaf, nfp, L_c, L_s
            
            !! input / output
            !integer :: n_ang
            !real(dp) :: find_f_plot(n_ang)
            !real(dp) :: dep_var(n_ang)
            !real(dp) :: par_var
            !real(dp) :: alpha
            
            !! local variables
            !integer :: jd
            !real(dp), allocatable :: cs(:,:,:)                                  ! (co)sines for all pol m and tor n
            !real(dp) :: lam
            
            !allocate(cs(0:mpol-1,-ntor:ntor,2))
            !do jd = 1, n_ang
                !if (use_tor_flux) then                                          ! looking for zeta (dep. var)
                    !! cosines and sines
                    !ierr = calc_ang_grid(cs,mpol,ntor,nfp,par_var,dep_var(jd))
                    !CHCKERR('')
                    !! lambda
                    !lam = f2r(L_c(:,:,kd,0),L_s(:,:,kd,0),cs,mpol,ntor,nfp,&
                        !&ierr=ierr)
                    !CHCKERR('')
                    !find_f_plot(jd) = dep_var(jd)-(par_var+lam)/iotaf(kd)-alpha
                !else                                                            ! looking for theta (dep. var)
                    !! cosines and sines
                    !ierr = calc_ang_grid(cs,mpol,ntor,nfp,dep_var(jd),par_var)
                    !CHCKERR('')
                    !! lambda
                    !lam = f2r(L_c(:,:,kd,0),L_s(:,:,kd,0),cs,mpol,ntor,nfp,&
                        !&ierr=ierr)
                    !CHCKERR('')
                    !find_f_plot(jd) = par_var-(dep_var(jd)+lam)/iotaf(kd)-alpha
                !end if
            !end do
        !end function find_f_plot
    end function test_ang_B
    
    integer function test_add_arr_mult() result(ierr)
        use utilities, only: add_arr_mult
        use eq_ops, only: init_eq, calc_RZL, calc_flux_q
        use coord_ops, only: calc_eqd_grid
        use eq_vars, only: R_E, Z_E, theta => theta_E, &
            &zeta => zeta_E, q_saf_E, pres_E, grp_n_r_eq
        use X_vars, only: n_par
        !use num_vars, only: max_deriv
        
        character(*), parameter :: rout_name = 'test_add_arr_mult'
        
        ! local variables
        integer :: case_nr
        
        ! initialize ierr
        ierr = 0
        
        write(*,*) 'test add_arr_mult?'
        if(yes_no(.false.)) then
            output_i = 0
            write(*,*) 'Testing whether add_arr_mult correctly multiplies two &
                &functions, taking into account derivatives'
            do 
                case_nr = 1
                write(*,*) 'Choose which case to test'
                write(*,*) '    1. arr_1 in 3 coords and arr_2 in 3 coords [def]'
                write(*,*) '    2. arr_1 in 3 coords and arr_2 in 1 coords'
                write(*,*) '    3. arr_1 in 1 coords and arr_2 in 1 coords'
                read(*,*) case_nr
                if (case_nr.eq.1 .or. case_nr.eq.2 .or. case_nr.eq.3) then
                    exit
                else
                    write(*,*) 'You have to choose one of the three &
                        &possibilities...'
                end if
            end do
            
            select case (case_nr)
                case(1)
                    ierr = case_3_3()
                    CHCKERR('')
                case(2)
                    ierr = case_3_1()
                    CHCKERR('')
                case(3)
                    ierr = case_1_1()
                    CHCKERR('')
                case default
                    write(*,*) 'ERROR: no case number associated with '&
                        &//trim(i2str(case_nr))
                    write(*,*) 'How did you get here???'
                    stop
            end select
            
            
            write(*,*) 'Stopping'
            stop
        end if
    contains
        integer function case_3_3() result(ierr)
            character(*), parameter :: rout_name = 'case_3_3'
            
            ! local variables
            real(dp) :: RZ(n_par,grp_n_r_eq), RZ_num(grp_n_r_eq)
            integer :: id, jd, kd
            
            ! initialize ierr
            ierr = 0
            
            ! calculate equilibrium quantities
            write(*,*) 'calculate RZL'
            if (allocated(zeta)) deallocate(zeta)
            allocate(zeta(n_par,grp_n_r_eq)); zeta = 0.0_dp
            if (allocated(theta)) deallocate(theta)
            allocate(theta(n_par,grp_n_r_eq)); theta = 0.0_dp
            
            ierr = init_eq()
            
            zeta = 0.4*pi/2
            do jd = 1,grp_n_r_eq
                ierr = calc_eqd_grid(theta(:,jd),n_par,0.0_dp*pi,3.0_dp*pi)
                CHCKERR('')
            end do
            
            do kd = 0, 3
                do jd = 0, 3
                    do id = 0, 3
                        ierr = calc_RZL([id,jd,kd])
                        CHCKERR('')
                    end do
                end do
            end do
            
            ! multiply
            write(*,*) 'multiply R and Z'
            RZ = 0.0_dp
            ierr = add_arr_mult(R_E,Z_E,RZ,[0,0,0])
            CHCKERR('')
            RZ_num = R_E(5,:,0,0,0)*Z_E(5,:,0,0,0)
            call print_GP_2D('RZ (calc,num) at par = 5','',&
                &reshape([RZ(5,:),RZ_num],[grp_n_r_eq,2]))
            call print_GP_2D('diff RZ (calc-num) at par = 5','',&
                &RZ(5,:)-RZ_num)
            
            ! derive multiplied values
            write(*,*) 'derive RZ'
            ! Dr
            RZ = 0.0_dp
            ierr = add_arr_mult(R_E,Z_E,RZ,[1,0,0])
            CHCKERR('')
            RZ_num = R_E(5,:,1,0,0)*Z_E(5,:,0,0,0) + &
                &R_E(5,:,0,0,0)*Z_E(5,:,1,0,0)
            call print_GP_2D('DrRZ (calc,num) at par = 5','',&
                &reshape([RZ(5,:),RZ_num],[grp_n_r_eq,2]))
            call print_GP_2D('diff DrRZ (calc-num) at par = 5','',&
                &RZ(5,:)-RZ_num)
            ! Dtheta
            RZ = 0.0_dp
            ierr = add_arr_mult(R_E,Z_E,RZ,[0,1,0])
            CHCKERR('')
            RZ_num = R_E(5,:,0,1,0)*Z_E(5,:,0,0,0) + &
                &R_E(5,:,0,0,0)*Z_E(5,:,0,1,0)
            call print_GP_2D('DtRZ (calc,num) at par = 5','',&
                &reshape([RZ(5,:),RZ_num],[grp_n_r_eq,2]))
            call print_GP_2D('diff DtRZ (calc-num) at par = 5','',&
                &RZ(5,:)-RZ_num)
            ! Dzeta
            RZ = 0.0_dp
            ierr = add_arr_mult(R_E,Z_E,RZ,[0,0,1])
            CHCKERR('')
            RZ_num = R_E(5,:,0,0,1)*Z_E(5,:,0,0,0) + &
                &R_E(5,:,0,0,0)*Z_E(5,:,0,0,1)
            call print_GP_2D('DzRZ (calc,num) at par = 5','',&
                &reshape([RZ(5,:),RZ_num],[grp_n_r_eq,2]))
            call print_GP_2D('diff DzRZ (calc-num) at par = 5','',&
                &RZ(5,:)-RZ_num)
            
            ! double derive multiplied values
            write(*,*) 'double derive RZ'
            ! Drr
            RZ = 0.0_dp
            ierr = add_arr_mult(R_E,Z_E,RZ,[2,0,0])
            CHCKERR('')
            RZ_num = R_E(5,:,2,0,0)*Z_E(5,:,0,0,0) + &
                &2.0_dp * R_E(5,:,1,0,0)*Z_E(5,:,1,0,0) + &
                &R_E(5,:,0,0,0)*Z_E(5,:,2,0,0)
            call print_GP_2D('DrrRZ (calc,num) at par = 5','',&
                &reshape([RZ(5,:),RZ_num],[grp_n_r_eq,2]))
            call print_GP_2D('diff DrrRZ (calc-num) at par = 5','',&
                &RZ(5,:)-RZ_num)
            ! Drtheta
            RZ = 0.0_dp
            ierr = add_arr_mult(R_E,Z_E,RZ,[1,1,0])
            CHCKERR('')
            RZ_num = R_E(5,:,1,1,0)*Z_E(5,:,0,0,0) + &
                &R_E(5,:,1,0,0)*Z_E(5,:,0,1,0) + &
                &R_E(5,:,0,1,0)*Z_E(5,:,1,0,0) + &
                &R_E(5,:,0,0,0)*Z_E(5,:,1,1,0)
            call print_GP_2D('DrtRZ (calc,num) at par = 5','',&
                &reshape([RZ(5,:),RZ_num],[grp_n_r_eq,2]))
            call print_GP_2D('diff DrtRZ (calc-num) at par = 5','',&
                &RZ(5,:)-RZ_num)
            ! Drzeta
            RZ = 0.0_dp
            ierr = add_arr_mult(R_E,Z_E,RZ,[1,0,1])
            CHCKERR('')
            RZ_num = R_E(5,:,1,0,1)*Z_E(5,:,0,0,0) + &
                &R_E(5,:,1,0,0)*Z_E(5,:,0,0,1) + &
                &R_E(5,:,0,0,1)*Z_E(5,:,1,0,0) + &
                &R_E(5,:,0,0,0)*Z_E(5,:,1,0,1)
            call print_GP_2D('DrzRZ (calc,num) at par = 5','',&
                &reshape([RZ(5,:),RZ_num],[grp_n_r_eq,2]))
            call print_GP_2D('diff DrzRZ (calc-num) at par = 5','',&
                &RZ(5,:)-RZ_num)
            ! Dthetatheta
            RZ = 0.0_dp
            ierr = add_arr_mult(R_E,Z_E,RZ,[0,2,0])
            CHCKERR('')
            RZ_num = R_E(5,:,0,2,0)*Z_E(5,:,0,0,0) + &
                &2.0_dp * R_E(5,:,0,1,0)*Z_E(5,:,0,1,0) + &
                &R_E(5,:,0,0,0)*Z_E(5,:,0,2,0)
            call print_GP_2D('DttRZ (calc,num) at par = 5','',&
                &reshape([RZ(5,:),RZ_num],[grp_n_r_eq,2]))
            call print_GP_2D('diff DttRZ (calc-num) at par = 5','',&
                &RZ(5,:)-RZ_num)
            ! Dtzeta
            RZ = 0.0_dp
            ierr = add_arr_mult(R_E,Z_E,RZ,[0,1,1])
            CHCKERR('')
            RZ_num = R_E(5,:,0,1,1)*Z_E(5,:,0,0,0) + &
                &R_E(5,:,0,1,0)*Z_E(5,:,0,0,1) + &
                &R_E(5,:,0,0,1)*Z_E(5,:,0,1,0) + &
                &R_E(5,:,0,0,0)*Z_E(5,:,0,1,1)
            call print_GP_2D('DtzRZ (calc,num) at par = 5','',&
                &reshape([RZ(5,:),RZ_num],[grp_n_r_eq,2]))
            call print_GP_2D('diff DtzRZ (calc-num) at par = 5','',&
                &RZ(5,:)-RZ_num)
            ! Dzetazeta
            RZ = 0.0_dp
            ierr = add_arr_mult(R_E,Z_E,RZ,[0,0,2])
            CHCKERR('')
            RZ_num = R_E(5,:,0,0,2)*Z_E(5,:,0,0,0) + &
                &2.0_dp * R_E(5,:,0,0,1)*Z_E(5,:,0,0,1) + &
                &R_E(5,:,0,0,0)*Z_E(5,:,0,0,2)
            call print_GP_2D('DzzRZ (calc,num) at par = 5','',&
                &reshape([RZ(5,:),RZ_num],[grp_n_r_eq,2]))
            call print_GP_2D('diff DzzRZ (calc-num) at par = 5','',&
                &RZ(5,:)-RZ_num)
            
            ! higher order derive multiplied values
            write(*,*) 'higher order derive RZ'
            RZ = 0.0_dp
            ierr = add_arr_mult(R_E,Z_E,RZ,[1,2,3])
            CHCKERR('')
            RZ_num = R_E(5,:,1,2,3)*Z_E(5,:,0,0,0) + &
                &3.0_dp * R_E(5,:,1,2,2)*Z_E(5,:,0,0,1) + &
                &3.0_dp * R_E(5,:,1,2,1)*Z_E(5,:,0,0,2) + &
                &1.0_dp * R_E(5,:,1,2,0)*Z_E(5,:,0,0,3) + &
                &2.0_dp * R_E(5,:,1,1,3)*Z_E(5,:,0,1,0) + &
                &6.0_dp * R_E(5,:,1,1,2)*Z_E(5,:,0,1,1) + &
                &6.0_dp * R_E(5,:,1,1,1)*Z_E(5,:,0,1,2) + &
                &2.0_dp * R_E(5,:,1,1,0)*Z_E(5,:,0,1,3) + &
                &1.0_dp * R_E(5,:,1,0,3)*Z_E(5,:,0,2,0) + &
                &3.0_dp * R_E(5,:,1,0,2)*Z_E(5,:,0,2,1) + &
                &3.0_dp * R_E(5,:,1,0,1)*Z_E(5,:,0,2,2) + &
                &1.0_dp * R_E(5,:,1,0,0)*Z_E(5,:,0,2,3) + &
                &1.0_dp * R_E(5,:,0,2,3)*Z_E(5,:,1,0,0) + &
                &3.0_dp * R_E(5,:,0,2,2)*Z_E(5,:,1,0,1) + &
                &3.0_dp * R_E(5,:,0,2,1)*Z_E(5,:,1,0,2) + &
                &1.0_dp * R_E(5,:,0,2,0)*Z_E(5,:,1,0,3) + &
                &2.0_dp * R_E(5,:,0,1,3)*Z_E(5,:,1,1,0) + &
                &6.0_dp * R_E(5,:,0,1,2)*Z_E(5,:,1,1,1) + &
                &6.0_dp * R_E(5,:,0,1,1)*Z_E(5,:,1,1,2) + &
                &2.0_dp * R_E(5,:,0,1,0)*Z_E(5,:,1,1,3) + &
                &1.0_dp * R_E(5,:,0,0,3)*Z_E(5,:,1,2,0) + &
                &3.0_dp * R_E(5,:,0,0,2)*Z_E(5,:,1,2,1) + &
                &3.0_dp * R_E(5,:,0,0,1)*Z_E(5,:,1,2,2) + &
                &1.0_dp * R_E(5,:,0,0,0)*Z_E(5,:,1,2,3)
            call print_GP_2D('DrttzzzRZ (calc,num) at par = 5','',&
                &reshape([RZ(5,:),RZ_num],[grp_n_r_eq,2]))
            call print_GP_2D('diff DrttzzzRZ (calc-num) at par = 5','',&
                &RZ(5,:)-RZ_num)
        end function case_3_3
        
        integer function case_3_1() result(ierr)
            character(*), parameter :: rout_name = 'case_3_1'
            
            ! local variables
            real(dp) :: Rq(n_par,grp_n_r_eq), Rq_num(grp_n_r_eq)
            integer :: id, jd, kd
            
            ! initialize ierr
            ierr = 0
            
            ! calculate equilibrium quantities
            write(*,*) 'Calculate RZL'
            if (allocated(zeta)) deallocate(zeta)
            allocate(zeta(n_par,grp_n_r_eq)); zeta = 0.0_dp
            if (allocated(theta)) deallocate(theta)
            allocate(theta(n_par,grp_n_r_eq)); theta = 0.0_dp
            
            ierr = init_eq()
            
            zeta = 0.4*pi/2
            do jd = 1,grp_n_r_eq
                ierr = calc_eqd_grid(theta(:,jd),n_par, 0.0_dp*pi, 3.0_dp*pi)
            CHCKERR('')
            end do
            
            do kd = 0, 3
                do jd = 0, 3
                    do id = 0, 3
                        ierr = calc_RZL([id,jd,kd])
                        CHCKERR('')
                    end do
                end do
            end do
            
            write(*,*) 'calculate q_saf_E'
            ierr = calc_flux_q()
            CHCKERR('')
            
            ! multiply
            write(*,*) 'multiply R and q_saf_E'
            Rq = 0.0_dp
            ierr = add_arr_mult(R_E,q_saf_E,Rq,[0,0,0])
            CHCKERR('')
            Rq_num = R_E(5,:,0,0,0)*q_saf_E(:,0)
            call print_GP_2D('Rq (calc,num) at par = 5','',&
                &reshape([Rq(5,:),Rq_num],[grp_n_r_eq,2]))
            call print_GP_2D('diff Rq (calc-num) at par = 5','',&
                &Rq(5,:)-Rq_num)
            
            ! derive multiplied values
            write(*,*) 'derive Rq'
            ! Dr
            Rq = 0.0_dp
            ierr = add_arr_mult(R_E,q_saf_E,Rq,[1,0,0])
            CHCKERR('')
            Rq_num = R_E(5,:,1,0,0)*q_saf_E(:,0) + &
                &R_E(5,:,0,0,0)*q_saf_E(:,1)
            call print_GP_2D('DrRq (calc,num) at par = 5','',&
                &reshape([Rq(5,:),Rq_num],[grp_n_r_eq,2]))
            call print_GP_2D('diff DrRq (calc-num) at par = 5','',&
                &Rq(5,:)-Rq_num)
            ! Dtheta
            Rq = 0.0_dp
            ierr = add_arr_mult(R_E,q_saf_E,Rq,[0,1,0])
            CHCKERR('')
            Rq_num = R_E(5,:,0,1,0)*q_saf_E(:,0) + &
                &R_E(5,:,0,0,0)*0.0_dp
            call print_GP_2D('DtRq (calc,num) at par = 5','',&
                &reshape([Rq(5,:),Rq_num],[grp_n_r_eq,2]))
            call print_GP_2D('diff DtRq (calc-num) at par = 5','',&
                &Rq(5,:)-Rq_num)
            ! Dzeta
            Rq = 0.0_dp
            ierr = add_arr_mult(R_E,q_saf_E,Rq,[0,0,1])
            CHCKERR('')
            Rq_num = R_E(5,:,0,0,1)*q_saf_E(:,0) + &
                &R_E(5,:,0,0,0)*0.0_dp
            call print_GP_2D('DzRq (calc,num) at par = 5','',&
                &reshape([Rq(5,:),Rq_num],[grp_n_r_eq,2]))
            call print_GP_2D('diff DzRq (calc-num) at par = 5','',&
                &Rq(5,:)-Rq_num)
            
            ! double derive multiplied values
            write(*,*) 'double derive Rq'
            ! Drr
            Rq = 0.0_dp
            ierr = add_arr_mult(R_E,q_saf_E,Rq,[2,0,0])
            CHCKERR('')
            Rq_num = R_E(5,:,2,0,0)*q_saf_E(:,0) + &
                &2.0_dp * R_E(5,:,1,0,0)*q_saf_E(:,1) + &
                &R_E(5,:,0,0,0)*q_saf_E(:,2)
            call print_GP_2D('DrrRq (calc,num) at par = 5','',&
                &reshape([Rq(5,:),Rq_num],[grp_n_r_eq,2]))
            call print_GP_2D('diff DrrRq (calc-num) at par = 5','',&
                &Rq(5,:)-Rq_num)
            ! Drtheta
            Rq = 0.0_dp
            ierr = add_arr_mult(R_E,q_saf_E,Rq,[1,1,0])
            CHCKERR('')
            Rq_num = R_E(5,:,1,1,0)*q_saf_E(:,0) + &
                &R_E(5,:,1,0,0)*0.0_dp + &
                &R_E(5,:,0,1,0)*q_saf_E(:,1) + &
                &R_E(5,:,0,0,0)*0.0_dp
            call print_GP_2D('DrtRq (calc,num) at par = 5','',&
                &reshape([Rq(5,:),Rq_num],[grp_n_r_eq,2]))
            call print_GP_2D('diff DrtRq (calc-num) at par = 5','',&
                &Rq(5,:)-Rq_num)
            ! Drzeta
            Rq = 0.0_dp
            ierr = add_arr_mult(R_E,q_saf_E,Rq,[1,0,1])
            CHCKERR('')
            Rq_num = R_E(5,:,1,0,1)*q_saf_E(:,0) + &
                &R_E(5,:,1,0,0)*0.0_dp + &
                &R_E(5,:,0,0,1)*q_saf_E(:,1) + &
                &R_E(5,:,0,0,0)*0.0_dp
            call print_GP_2D('DrzRq (calc,num) at par = 5','',&
                &reshape([Rq(5,:),Rq_num],[grp_n_r_eq,2]))
            call print_GP_2D('diff DrzRq (calc-num) at par = 5','',&
                &Rq(5,:)-Rq_num)
            ! Dthetatheta
            Rq = 0.0_dp
            ierr = add_arr_mult(R_E,q_saf_E,Rq,[0,2,0])
            CHCKERR('')
            Rq_num = R_E(5,:,0,2,0)*q_saf_E(:,0) + &
                &2.0_dp * R_E(5,:,0,1,0)*0.0_dp + &
                &R_E(5,:,0,0,0)*0.0_dp
            call print_GP_2D('DttRq (calc,num) at par = 5','',&
                &reshape([Rq(5,:),Rq_num],[grp_n_r_eq,2]))
            call print_GP_2D('diff DttRq (calc-num) at par = 5','',&
                &Rq(5,:)-Rq_num)
            ! Dtzeta
            Rq = 0.0_dp
            ierr = add_arr_mult(R_E,q_saf_E,Rq,[0,1,1])
            CHCKERR('')
            Rq_num = R_E(5,:,0,1,1)*q_saf_E(:,0) + &
                &R_E(5,:,0,1,0)*0.0_dp + &
                &R_E(5,:,0,0,1)*0.0_dp + &
                &R_E(5,:,0,0,0)*0.0_dp
            call print_GP_2D('DtzRq (calc,num) at par = 5','',&
                &reshape([Rq(5,:),Rq_num],[grp_n_r_eq,2]))
            call print_GP_2D('diff DtzRq (calc-num) at par = 5','',&
                &Rq(5,:)-Rq_num)
            ! Dzetazeta
            Rq = 0.0_dp
            ierr = add_arr_mult(R_E,q_saf_E,Rq,[0,0,2])
            CHCKERR('')
            Rq_num = R_E(5,:,0,0,2)*q_saf_E(:,0) + &
                &2.0_dp * R_E(5,:,0,0,1)*0.0_dp + &
                &R_E(5,:,0,0,0)*0.0_dp
            call print_GP_2D('DzzRq (calc,num) at par = 5','',&
                &reshape([Rq(5,:),Rq_num],[grp_n_r_eq,2]))
            call print_GP_2D('diff DzzRq (calc-num) at par = 5','',&
                &Rq(5,:)-Rq_num)
            
            ! higher order derive multiplied values
            write(*,*) 'higher order derive Rq'
            Rq = 0.0_dp
            ierr = add_arr_mult(R_E,q_saf_E,Rq,[1,2,3])
            CHCKERR('')
            Rq_num = R_E(5,:,1,2,3)*q_saf_E(:,0) + &
                &3.0_dp * R_E(5,:,1,2,2)*0.0_dp + &
                &3.0_dp * R_E(5,:,1,2,1)*0.0_dp + &
                &1.0_dp * R_E(5,:,1,2,0)*0.0_dp + &
                &2.0_dp * R_E(5,:,1,1,3)*0.0_dp + &
                &6.0_dp * R_E(5,:,1,1,2)*0.0_dp + &
                &6.0_dp * R_E(5,:,1,1,1)*0.0_dp + &
                &2.0_dp * R_E(5,:,1,1,0)*0.0_dp + &
                &1.0_dp * R_E(5,:,1,0,3)*0.0_dp + &
                &3.0_dp * R_E(5,:,1,0,2)*0.0_dp + &
                &3.0_dp * R_E(5,:,1,0,1)*0.0_dp + &
                &1.0_dp * R_E(5,:,1,0,0)*0.0_dp + &
                &1.0_dp * R_E(5,:,0,2,3)*q_saf_E(:,1) + &
                &3.0_dp * R_E(5,:,0,2,2)*0.0_dp + &
                &3.0_dp * R_E(5,:,0,2,1)*0.0_dp + &
                &1.0_dp * R_E(5,:,0,2,0)*0.0_dp + &
                &2.0_dp * R_E(5,:,0,1,3)*0.0_dp + &
                &6.0_dp * R_E(5,:,0,1,2)*0.0_dp + &
                &6.0_dp * R_E(5,:,0,1,1)*0.0_dp + &
                &2.0_dp * R_E(5,:,0,1,0)*0.0_dp + &
                &1.0_dp * R_E(5,:,0,0,3)*0.0_dp + &
                &3.0_dp * R_E(5,:,0,0,2)*0.0_dp + &
                &3.0_dp * R_E(5,:,0,0,1)*0.0_dp + &
                &1.0_dp * R_E(5,:,0,0,0)*0.0_dp
            call print_GP_2D('DrttzzzRq (calc,num) at par = 5','',&
                &reshape([Rq(5,:),Rq_num],[grp_n_r_eq,2]))
            call print_GP_2D('diff DrttzzzRq (calc-num) at par = 5','',&
                &Rq(5,:)-Rq_num)
        end function case_3_1
        
        integer function case_1_1() result(ierr)
            character(*), parameter :: rout_name = 'case_1_1'
            
            ! local variables
            real(dp) :: pq(grp_n_r_eq), pq_num(grp_n_r_eq)
            integer :: jd
            
            ! initialize ierr
            ierr = 0
            
            ! calculate equilibrium quantities
            write(*,*) 'calculate RZL'
            if (allocated(zeta)) deallocate(zeta)
            allocate(zeta(n_par,grp_n_r_eq)); zeta = 0.0_dp
            if (allocated(theta)) deallocate(theta)
            allocate(theta(n_par,grp_n_r_eq)); theta = 0.0_dp
            
            ierr = init_eq()
            
            zeta = 0.4*pi/2
            do jd = 1,grp_n_r_eq
                ierr = calc_eqd_grid(theta(:,jd),n_par,0.0_dp*pi,3.0_dp*pi)
                CHCKERR('')
            end do
            
            write(*,*) 'calculate q_saf_E and pres_E'
            ierr = calc_flux_q()
            CHCKERR('')
            
            ! multiply
            write(*,*) 'multiply pres_E and q_saf_E'
            pq = 0.0_dp
            ierr = add_arr_mult(pres_E,q_saf_E,pq,[0,0,0])
            CHCKERR('')
            pq_num = pres_E(:,0)*q_saf_E(:,0)
            call print_GP_2D('pq (calc,num) = 5','',&
                &reshape([pq,pq_num],[grp_n_r_eq,2]))
            call print_GP_2D('diff pq (calc-num) = 5','',&
                &pq-pq_num)
            
            ! derive multiplied values
            write(*,*) 'derive pq'
            ! Dr
            pq = 0.0_dp
            ierr = add_arr_mult(pres_E,q_saf_E,pq,[1,0,0])
            CHCKERR('')
            pq_num = pres_E(:,1)*q_saf_E(:,0) + &
                &pres_E(:,0)*q_saf_E(:,1)
            call print_GP_2D('Drpq (calc,num) = 5','',&
                &reshape([pq,pq_num],[grp_n_r_eq,2]))
            call print_GP_2D('diff Drpq (calc-num) = 5','',&
                &pq-pq_num)
            ! Dtheta
            pq = 0.0_dp
            ierr = add_arr_mult(pres_E,q_saf_E,pq,[0,1,0])
            CHCKERR('')
            pq_num = 0.0_dp*q_saf_E(:,0) + &
                &pres_E(:,0)*0.0_dp
            call print_GP_2D('Dtpq (calc,num) = 5','',&
                &reshape([pq,pq_num],[grp_n_r_eq,2]))
            call print_GP_2D('diff Dtpq (calc-num) = 5','',&
                &pq-pq_num)
            ! Dzeta
            pq = 0.0_dp
            ierr = add_arr_mult(pres_E,q_saf_E,pq,[0,0,1])
            CHCKERR('')
            pq_num = 0.0_dp*q_saf_E(:,0) + &
                &pres_E(:,0)*0.0_dp
            call print_GP_2D('Dzpq (calc,num) = 5','',&
                &reshape([pq,pq_num],[grp_n_r_eq,2]))
            call print_GP_2D('diff Dzpq (calc-num) = 5','',&
                &pq-pq_num)
            
            ! double derive multiplied values
            write(*,*) 'double derive pq'
            ! Drr
            pq = 0.0_dp
            ierr = add_arr_mult(pres_E,q_saf_E,pq,[2,0,0])
            CHCKERR('')
            pq_num = pres_E(:,2)*q_saf_E(:,0) + &
                &2.0_dp * pres_E(:,1)*q_saf_E(:,1) + &
                &pres_E(:,0)*q_saf_E(:,2)
            call print_GP_2D('Drrpq (calc,num) = 5','',&
                &reshape([pq,pq_num],[grp_n_r_eq,2]))
            call print_GP_2D('diff Drrpq (calc-num) = 5','',&
                &pq-pq_num)
            ! Drtheta
            pq = 0.0_dp
            ierr = add_arr_mult(pres_E,q_saf_E,pq,[1,1,0])
            CHCKERR('')
            pq_num = 0.0_dp
            call print_GP_2D('Drtpq (calc,num) = 5','',&
                &reshape([pq,pq_num],[grp_n_r_eq,2]))
            call print_GP_2D('diff Drtpq (calc-num) = 5','',&
                &pq-pq_num)
            ! Drzeta
            pq = 0.0_dp
            ierr = add_arr_mult(pres_E,q_saf_E,pq,[2,0,1])
            CHCKERR('')
            pq_num = 0.0_dp
            call print_GP_2D('Drzpq (calc,num) = 5','',&
                &reshape([pq,pq_num],[grp_n_r_eq,2]))
            call print_GP_2D('diff Drzpq (calc-num) = 5','',&
                &pq-pq_num)
            ! Dthetatheta
            pq = 0.0_dp
            ierr = add_arr_mult(pres_E,q_saf_E,pq,[0,2,0])
            CHCKERR('')
            pq_num = 0.0_dp
            call print_GP_2D('Dttpq (calc,num) = 5','',&
                &reshape([pq,pq_num],[grp_n_r_eq,2]))
            call print_GP_2D('diff Dttpq (calc-num) = 5','',&
                &pq-pq_num)
            ! Dtzeta
            pq = 0.0_dp
            ierr = add_arr_mult(pres_E,q_saf_E,pq,[0,1,1])
            CHCKERR('')
            pq_num = 0.0_dp
            call print_GP_2D('Dtzpq (calc,num) = 5','',&
                &reshape([pq,pq_num],[grp_n_r_eq,2]))
            call print_GP_2D('diff Dtzpq (calc-num) = 5','',&
                &pq-pq_num)
            ! Dzetazeta
            pq = 0.0_dp
            ierr = add_arr_mult(pres_E,q_saf_E,pq,[0,0,2])
            CHCKERR('')
            pq_num = 0.0_dp
            call print_GP_2D('Dzzpq (calc,num) = 5','',&
                &reshape([pq,pq_num],[grp_n_r_eq,2]))
            call print_GP_2D('diff Dzzpq (calc-num) = 5','',&
                &pq-pq_num)
            
            ! higher order derive multiplied values
            write(*,*) 'higher order derive pq'
            pq = 0.0_dp
            ierr = add_arr_mult(pres_E,q_saf_E,pq,[1,2,3])
            CHCKERR('')
            pq_num = 0.0_dp
            call print_GP_2D('Drttzzzpq (calc,num) = 5','',&
                &reshape([pq,pq_num],[grp_n_r_eq,2]))
            call print_GP_2D('diff Drttzzzpq (calc-num) = 5','',&
                &pq-pq_num)
        end function case_1_1
    end function test_add_arr_mult
    
    ! test SLEPC
    ! (from $PETSC_DIR/src/eps/examples/tutorials/ex1f90.F90)
    integer function test_slepc() result(ierr)
#include <finclude/slepcepsdef.h>
        use slepceps
        
        character(*), parameter :: rout_name = 'test_slepc'
        
        ! local variables
        Mat :: A
        EPS :: solver
        PetscMPIInt :: rank, n_procs
        PetscInt :: n_m_X, n_r_X, n_block
        PetscInt :: i, j, Istart, Iend
        PetscInt :: nev, ncv, mpd
        EPSType :: tname
        PetscScalar, allocatable :: loc_block(:,:)
        PetscInt, allocatable :: loc_x(:), loc_y(:)
        
        ! initialize ierr
        ierr = 0
        
        write(*,*) 'test slepc?'
        if(yes_no(.false.)) then
            
            call SlepcInitialize(PETSC_NULL_CHARACTER,ierr)
            CHCKERR('')
            call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr)
            CHCKERR('')
            
            ! set options
            call PetscOptionsSetValue('-draw_pause','-1',ierr)
            CHCKERR('')
            call PetscOptionsSetValue('-eps_view','',ierr)
            CHCKERR('')
            
#if defined(PETSC_USE_COMPLEX)
            call PetscPrintf(PETSC_COMM_WORLD,"using complex numbers\n",ierr);
            CHCKERR('')
#else
            call PetscPrintf(PETSC_COMM_WORLD,"can only use complex numbers!\n",ierr);
            CHCKERR('')
            call SlepcFinalize(ierr)
            CHCKERR('')
            stop
#endif
            
            ! divide the n_r_X normal points between the n_procs processors
            n_m_X = 6                                                           ! size of block (=M)
            n_r_X = 10                                                          ! number of radial points
            call MPI_Comm_size(PETSC_COMM_WORLD,n_procs,ierr)                   ! number of processors
            CHCKERR('')
            if (n_procs.gt.n_r_X) then
                call PetscPrintf(PETSC_COMM_WORLD,"n_procs too high!\n",ierr);
                CHCKERR('')
                call SlepcFinalize(ierr)
                CHCKERR('')
                stop
            end if
            n_block = n_r_X/n_procs                                             ! number of radial points on this processor
            if (mod(n_r_X,n_procs).gt.0) then
                if (mod(n_r_X,n_procs).gt.rank) n_block = n_block + 1           ! add a point to the first ranks if there is a remainder
            end if
            
            ! create  a  matrix A  with the  appropriate number  of preallocated
            ! entries
            call MatCreateAIJ(PETSC_COMM_WORLD,n_block*n_m_X,n_block*n_m_X,&
                &n_r_X*n_m_X,n_r_X*n_m_X,3*n_m_X,PETSC_NULL_INTEGER,&
                &3*n_m_X,PETSC_NULL_INTEGER,A,ierr)
            CHCKERR('')
            
            ! fill the matrix A
            call MatGetOwnershipRange(A,Istart,Iend,ierr)
            CHCKERR('')
            !write(*,*) 'rank, Istart, Iend = ', rank, Istart, Iend
            allocate(loc_block(n_m_X,n_m_X))
            allocate(loc_x(n_m_X))
            allocate(loc_y(n_m_X))
            do i = 1,n_block
                if (rank.ne.0 .or. i.ne.1) then                                 ! don't do the left block for first normal point
                    ! for each to the left block add -1
                    loc_block = -1.0
                    loc_x = [(j, j = Istart,Iend-1)] + (i-1)*n_m_X
                    loc_y = loc_x-n_m_X
                    call MatSetValues(A,n_m_X,loc_x,n_m_X,loc_y,loc_block,&
                        &INSERT_VALUES,ierr)
                    CHCKERR('')
                end if
                ! for each block add imaginary unit PETSC_i
                loc_block = PETSC_i
                loc_x = [(j, j = Istart,Iend-1)] + (i-1)*n_m_X
                loc_y = loc_x
                call MatSetValues(A,n_m_X,loc_x,n_m_X,loc_y,loc_block,&
                    &INSERT_VALUES,ierr)
                CHCKERR('')
                if (rank.ne.n_procs-1 .or. i.ne.n_block) then                   ! don't do the right block for last normal point
                    ! for each block add 1
                    loc_block = +1.0
                    loc_x = [(j, j = Istart,Iend-1)] + (i-1)*n_m_X
                    loc_y = loc_x+n_m_X
                    call MatSetValues(A,n_m_X,loc_x,n_m_X,loc_y,loc_block,&
                        &INSERT_VALUES,ierr)
                    CHCKERR('')
                end if
            end do
            
            ! assemble the matrix and view it
            call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
            CHCKERR('')
            call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
            CHCKERR('')
            
            call MatView(A,PETSC_VIEWER_STDOUT_WORLD,ierr)
            CHCKERR('')
            call MatView(A,PETSC_VIEWER_DRAW_WORLD,ierr)
            CHCKERR('')

            ! create vector
            !call MatGetVecs(A,X,Y,ierr)
            CHCKERR('')
            !call VecGetOwnershipRange(X,Istart,Iend,ierr)
            CHCKERR('')
            !do i= Istart,Iend-1
                !val = rank*1000 + 100.0*i + 10.0*PETSC_i
                !call VecSetValues(X,one,i,val,ADD_VALUES,ierr)
                CHCKERR('')
            !end do
            !call VecAssemblyBegin(X,ierr)
            CHCKERR('')
            !call VecAssemblyEnd(X,ierr)
            CHCKERR('')
            
            ! matmult with matrix A
            !call VecDuplicate(X,Y,ierr)
            CHCKERR('')
            !call MatMult(A,X,Y,ierr)
            CHCKERR('')
              
            !call PetscObjectSetName(X, "starting vector X:",ierr)
            CHCKERR('')
            !call PetscObjectSetName(Y, "multiplied vector Y:",ierr)
            CHCKERR('')
            !call VecView(X,PETSC_VIEWER_STDOUT_WORLD,ierr)
            CHCKERR('')
            !call VecView(Y,PETSC_VIEWER_STDOUT_WORLD,ierr)
            CHCKERR('')
            !call VecDestroy(X,ierr)
            CHCKERR('')
            !call VecDestroy(Y,ierr)
            CHCKERR('')
            
            ! solve an EV problem
            call EPSCreate(PETSC_COMM_WORLD,solver,ierr)
            CHCKERR('')
            call EPSSetOperators(solver,A,PETSC_NULL_OBJECT,ierr)
            CHCKERR('')
            call EPSSetProblemType(solver,EPS_HEP,ierr)
            CHCKERR('')
            call EPSSetFromOptions(solver,ierr)
            CHCKERR('')
            call EPSSolve(solver,ierr) 
            CHCKERR('')

            ! get some information from the solver and display it
            call EPSGetType(solver,tname,ierr)
            CHCKERR('')
            if (rank .eq. 0) write(*,'(" solution method: ",A)') tname
            call EPSGetDimensions(solver,nev,ncv,mpd,ierr)
            CHCKERR('')
            if (rank .eq. 0) write(*,'(" number of req. EV: ",I4)') nev
            if (rank .eq. 0) write(*,'(" max. dim. of subspace: ",I4)') ncv
            if (rank .eq. 0) write(*,'(" max. dim. for projected problem: ",&
                &I4)') mpd

            ! destroy and finalize
            call EPSPrintSolution(solver,PETSC_NULL_OBJECT,ierr)
            CHCKERR('')
            call EPSDestroy(solver,ierr)
            CHCKERR('')
            call MatDestroy(A,ierr)
            CHCKERR('')
            
            call SlepcFinalize(ierr)
            CHCKERR('')
            
            
            write(*,*) 'Stopping'
            stop
        end if
    end function test_slepc
end module test
