!------------------------------------------------------------------------------!
!   Calculate the matrix  elements due to the Plasma of  the system of equations
!   that has to be solved as described in [ADD REF]
!------------------------------------------------------------------------------!
module X_ops
#include <PB3D_macros.h>
    use num_vars, only: dp, iu, max_str_ln, pi
    use output_ops, only: lvl_ud, writo, print_ar_2, print_GP_2D
    use str_ops, only: i2str, r2strt

    implicit none
    private
    public prepare_X, solve_EV_system, plot_X_vec

contains
    ! prepare the matrix elements by calculating KV and PV, which then will have
    ! to be integrated, with a complex exponential weighting function
    subroutine prepare_X
        use X_vars, only: init_X, calc_PV, calc_KV, calc_U, calc_extra, &
            &calc_V_int, dealloc_X_vars, &
            &PV0, PV1, PV2, KV0, KV1, KV2, PV_int, KV_int
        
        call writo('Calculating table of field line averages &
            &<V e^i(k-m)theta>')
        
        call lvl_ud(1)
        
        ! initialize the variables
        call writo('Initalizing variables...')
        call lvl_ud(1)
        call init_X
        call lvl_ud(-1)
        
        ! calculate U and DU
        call writo('Calculating U and DU...')
        call lvl_ud(1)
        call calc_U
        call lvl_ud(-1)
        
        ! calculate extra equilibrium quantities
        call writo('Calculating extra equilibrium quantities...')
        call lvl_ud(1)
        call calc_extra
        call lvl_ud(-1)
        
        ! Calculate PV0, PV1  and PV2 for all (k,m) pairs  and n_r (equilibrium)
        ! values of the normal coordinate
        call writo('Calculating PV0, PV1 and PV2...')
        call lvl_ud(1)
        call calc_PV
        call lvl_ud(-1)
        
        ! Calculate KV0, KV1  and KV2 for all (k,m) pairs  and n_r (equilibrium)
        ! values of the normal coordinate
        call writo('Calculating KV0, KV1 and KV2...')
        call lvl_ud(1)
        call calc_KV
        call lvl_ud(-1)
        
        ! Calculate PV_int = <PVi e^(k-m)theta>
        call writo('Taking field average of PV')
        call lvl_ud(1)
        call calc_V_int(PV0,PV_int(:,:,:,1))
        call calc_V_int(PV1,PV_int(:,:,:,2))
        call calc_V_int(PV2,PV_int(:,:,:,3))
        call lvl_ud(-1)
        
        ! Calculate KV_int = <KVi e^(k-m)theta>
        call writo('Taking field average of KV')
        call lvl_ud(1)
        call calc_V_int(KV0,KV_int(:,:,:,1))
        call calc_V_int(KV1,KV_int(:,:,:,2))
        call calc_V_int(KV2,KV_int(:,:,:,3))
        call lvl_ud(-1)
        
        ! deallocate equilibrium variables
        call writo('deallocating unused variables')
        call dealloc_X_vars
        
        call lvl_ud(-1)
        
        call writo('Done calculating')
    end subroutine prepare_X
    
    ! set-up and  solve the  EV system  by discretizing  the equations  in n_r_X
    ! normal points,  making use of  PV0, PV1 and  PV2, interpolated in  the n_r
    ! (equilibrium) values
    integer function solve_EV_system(A_terms,B_terms,A_info,B_info) result(ierr)
        use num_vars, only: EV_style
        use str_ops, only: i2str
        use slepc_ops, only: solve_EV_system_slepc
        
        character(*), parameter :: rout_name = 'solve_EV_system'
        
        ! input / output
        complex(dp), intent(inout), allocatable, optional :: A_terms(:,:,:,:)   ! termj_int of matrix A from previous Richardson loop
        complex(dp), intent(inout), allocatable, optional :: B_terms(:,:,:,:)   ! termj_int of matrix B from previous Richardson loop
        real(dp), intent(inout), optional :: A_info(2)                          ! info about A_terms: min of r and step_size
        real(dp), intent(inout), optional :: B_info(2)                          ! info about B_terms: min of r and step_size
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        select case (EV_style)
            case(1)                                                             ! slepc solver for EV problem
                ! solve the system
                ierr = solve_EV_system_slepc(A_terms,B_terms,A_info,B_info)
                CHCKERR('')
                
            case default
                err_msg = 'No EV solver style associated with '//&
                    &trim(i2str(EV_style))
                ierr = 1
                CHCKERR(err_msg)
        end select
    end function solve_EV_system
    
    ! plots an eigenfunction in either 1D (normal coord.) or 2D (normal coord. &
    ! par coord.), depending on whether the variable theta is provided
    ! The  individual modes are  plotted at t  = 0 and  t = (pi/2)/omega  and an
    ! animated plot is produced as well in a .gif output file
    ! [MPI] Collective call
    integer function plot_X_vec(X_vec,X_val,X_id,theta) result(ierr)
        use num_vars, only: grp_rank, min_r_X, max_r_X, max_n_plots
        use output_ops, only: draw_GP, draw_GP_animated
        use X_vars, only: n_r_X, n_m_X, n_X, m_X
        use MPI_ops, only: get_ser_X_vec
        use eq_vars, only: q_saf_full
        
        character(*), parameter :: rout_name = 'plot_X_vec'
        
        ! input / output
        complex(dp), intent(in) :: X_vec(n_m_X,n_r_X)                           ! MPI Eigenvector
        complex(dp), intent(in) :: X_val                                        ! Eigenvalue
        integer, intent(in) :: X_id                                             ! nr. of Eigenvalue
        real(dp), intent(in), optional :: theta(:)                              ! the parallel range [2pi]
        
        ! local variables
        real(dp), allocatable :: x_plot(:,:,:)                                  ! x_axis of plot
        real(dp), allocatable :: y_plot(:,:,:)                                  ! y values of plot
        complex(dp), allocatable :: ser_X_vec(:,:)                              ! serial version of X_vec
        integer :: id, jd, kd                                                   ! counter
        real(dp), allocatable :: max_of_modes(:)                                ! maximum of each mode
        real(dp), allocatable :: max_of_modes_r(:)                              ! flux surface where max of mode occurs
        complex(dp) :: omega                                                    ! eigenvalue
        integer :: n_t                                                          ! nr. of time steps in quarter period
        integer :: n_n_t                                                        ! how many quarter periods
        character(len=max_str_ln) :: plot_title                                 ! title for plots
        real(dp) :: ranges(3,2)                                                 ! the range of the plots
        
        ! initialize ierr
        ierr = 0
        
        ! initialize other quantities
        if (grp_rank.eq.0) then                                                 ! only for group masters
            allocate(ser_X_vec(n_m_X,n_r_X))
        else
            allocate(ser_X_vec(0,0))
        end if
        
        ! get serial version of X_vec on group master
        ierr = get_ser_X_vec(X_vec,ser_X_vec,n_m_X)
        CHCKERR('')
        
        call lvl_ud(1)
        
        ! calculations for group master
        if (grp_rank.eq.0) then                                                 ! only for group masters
            ! get omega
            omega = sqrt(X_val)
            
            ! if  omega  is  predominantly  real,  the  Eigenfunction  explodes,
            ! so  limit  n_n_t to  1.  If  it  is predominantly  imaginary,  the
            ! Eigenfunction oscillates, so choose n_n_t = 8 for 2 whole periods
            if (abs(realpart(omega)/imagpart(omega)).gt.1) then
                n_n_t = 1                                                       ! 1 quarter period
                n_t = 20                                                        ! 20 points per quarter period
            else
                n_n_t = 8                                                       ! 8 quarter periods
                n_t = 5                                                         ! 5 points per quarter period
            end if
            
            ! set up x-axis values
            allocate(x_plot(1:n_r_X,1:n_m_X,n_n_t*n_t))
            do id = 1,n_r_X
                x_plot(id,:,:) = min_r_X + (id-1.0)/(n_r_X-1)*(max_r_X-min_r_X)
            end do
            
            ! set up y-axis values (not same for all time steps)
            ! allocate variables
            allocate(y_plot(n_r_X,n_m_X,n_n_t*n_t))
            allocate(max_of_modes(n_m_X))
            allocate(max_of_modes_r(n_m_X))
            max_of_modes = 0.0_dp
            max_of_modes_r = 0.0_dp
            
            ! check whether 2D or 3D plot
            if (present(theta)) then                                            ! 3D plot
            else                                                                ! 2D plot
                ! plot for all time steps the Eigenfunction at theta = 0
                do id = 1,n_n_t*n_t
                    ! loop over all normal points of all modes
                    do kd = 1,n_r_X
                        ! loop over all modes
                        do jd = 1,n_m_X
                            ! check for maximum of mode jd and normal point jd
                            y_plot(kd,jd,id) = realpart(ser_X_vec(jd,kd) * &
                                &exp(omega/abs(omega)*(id-1)/(n_t-1)*2*pi))
                            if (id.eq.1 .or. id.eq.n_t+1) then
                                if (y_plot(kd,jd,id).gt.max_of_modes(jd)) then
                                    max_of_modes(jd) = y_plot(kd,jd,id)
                                    max_of_modes_r(jd) = min_r_X + &
                                        &(max_r_X-min_r_X)*(kd-1.0)/(n_r_X-1.0)
                                end if
                            end if
                        end do
                    end do
                    if (id.eq.1 .or. id.eq.n_t+1) then
                        !! plot places of maximum
                        if (id.eq.1) then
                            !call writo('Information about the maximum of the &
                                !&modes for time t=0')
                            plot_title = 'Eigenvalue '//trim(i2str(X_id))//&
                                &' - Eigenvector (t = 0)'
                        else if (id.eq.n_t+1) then
                            !call writo('Information about the maximum of the &
                                !&modes for time t=T_0/4')
                            plot_title = 'Eigenvalue '//trim(i2str(X_id))//&
                                &' - Eigenvector (t = T_0/4)'
                        end if
                        !call print_GP_2D('maximum of the modes','',max_of_modes)
                        !call print_GP_2D('place of maximum of the modes','',&
                            !&max_of_modes_r)
                        ! plot the modes on screen
                        call print_GP_2D(plot_title,'Eigenvector.dat',&
                            &y_plot(:,:,id),x=x_plot(:,:,id),draw=.true.)
                        ! plot them in file as well
                        call draw_GP(plot_title,'Eigenvector.dat',&
                            &n_m_X,.true.,.false.)
                    end if
                end do
            end if
            
            deallocate(max_of_modes)
            deallocate(max_of_modes_r)
            
            if (n_m_X.gt.1) then
                ! plot the individual modes in time
                do jd = 1,min(n_m_X,max_n_plots)
                    ! set plot_title and ranges
                    plot_title = 'Eigenvalue '//trim(i2str(X_id))//' - &
                        &mode '//trim(i2str(jd))//' of '//trim(i2str(n_m_X))
                    ranges(1,:) = [minval(x_plot(:,jd,:)),&
                        &maxval(x_plot(:,jd,:))]
                    ranges(2,:) = [minval(y_plot(:,jd,:)),&
                        &maxval(y_plot(:,jd,:))]
                    call print_GP_2D(trim(plot_title),'mode'//trim(i2str(jd)),&
                        &y_plot(:,jd,:),x=x_plot(:,jd,:),draw=.false.)
                    call draw_GP_animated(trim(plot_title),'mode'//&
                        &trim(i2str(jd)),n_n_t*n_t,.true.,ranges(1:2,:))
                end do
                if (n_m_X.eq.max_n_plots+1) then
                    call writo('WARNING: Skipping plot for mode '//&
                        &trim(i2str(max_n_plots+1)))
                else if (n_m_X.gt.max_n_plots+1) then
                    call writo('WARNING: Skipping plots for modes '//&
                        &trim(i2str(max_n_plots+1))//'..'//trim(i2str(n_m_X)))
                end if
            end if
            
            ! plot the global mode in time
            do jd = 2,n_m_X
                y_plot(:,1,:) = y_plot(:,1,:) + y_plot(:,jd,:)
            end do
            plot_title = 'Eigenvalue '//trim(i2str(X_id))//' - global mode'
            ranges(1,:) = [minval(x_plot(:,1,:)),maxval(x_plot(:,1,:))]
            ranges(2,:) = [minval(y_plot(:,1,:)),maxval(y_plot(:,1,:))]
            call print_GP_2D(trim(plot_title),'global mode',&
                &y_plot(:,1,:),x=x_plot(:,1,:),draw=.false.)
            call draw_GP_animated(trim(plot_title),'global mode',&
                &n_n_t*n_t,.true.,ranges(1:2,:))
            
            ! deallocate variables
            deallocate(ser_X_vec)
            deallocate(x_plot)
        end if
        
        call lvl_ud(-1)
    end function plot_X_vec
end module
