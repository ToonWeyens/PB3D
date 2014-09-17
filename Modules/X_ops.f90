!------------------------------------------------------------------------------!
!   Calculate the matrix  elements due to the Plasma of  the system of equations
!   that has to be solved as described in [ADD REF]
!------------------------------------------------------------------------------!
module X_ops
#include <PB3D_macros.h>
    use num_vars, only: dp, iu, max_str_ln, pi
    use output_ops, only: lvl_ud, writo, print_ar_2, print_GP_2D, print_GP_3D
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
    integer function solve_EV_system(use_guess) result(ierr)
        use num_vars, only: EV_style
        use str_ops, only: i2str
        use slepc_ops, only: solve_EV_system_slepc
        
        character(*), parameter :: rout_name = 'solve_EV_system'
        
        ! input / output
        logical, intent(in) :: use_guess                                        ! whether to use a guess or not
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        select case (EV_style)
            case(1)                                                             ! slepc solver for EV problem
                ! solve the system
                ierr = solve_EV_system_slepc(use_guess)
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
        use num_vars, only: grp_rank, min_r_X, max_r_X, max_n_plots, &
            &alpha_job_nr
        use output_ops, only: draw_GP, draw_GP_animated
        use X_vars, only: n_r_X, n_m_X, n_X, m_X
        use VMEC_vars, only: n_r_eq
        use MPI_ops, only: get_ser_X_vec
        use eq_vars, only: q_saf_full
        use utilities, only: calc_interp, dis2con
        
        character(*), parameter :: rout_name = 'plot_X_vec'
        
        ! input / output
        complex(dp), intent(in) :: X_vec(n_m_X,n_r_X)                           ! MPI Eigenvector
        complex(dp), intent(in) :: X_val                                        ! Eigenvalue
        integer, intent(in) :: X_id                                             ! nr. of Eigenvalue (for output name)
        real(dp), intent(in), optional :: theta(2)                              ! the parallel range [2pi]
        
        ! local variables
        real(dp), allocatable :: x_plot(:,:,:,:)                                ! x_axis of plot
        real(dp), allocatable :: y_plot(:,:,:,:)                                ! y values of plot
        real(dp), allocatable :: z_plot(:,:,:,:)                                ! z values of plot
        complex(dp), allocatable :: ser_X_vec(:,:)                              ! serial version of X_vec
        integer :: id, jd, kd, ld                                               ! counter
        real(dp), allocatable :: max_of_modes(:)                                ! maximum of each mode
        real(dp), allocatable :: max_of_modes_x(:)                              ! flux surface where max of mode occurs
        real(dp), allocatable :: max_of_modes_y(:)                              ! parallel point where max of mode occurs
        complex(dp) :: omega                                                    ! eigenvalue
        integer :: n_t                                                          ! nr. of time steps in quarter period
        integer :: n_n_t                                                        ! how many quarter periods
        integer :: n_theta = 50                                                 ! how many parallel points
        integer :: max_n_norm = 50                                              ! max nr of normal points to plot
        integer :: n_norm                                                       ! nr. of normal points
        real(dp) :: inc_norm                                                    ! normal increment for plots
        real(dp) :: kd_loc                                                      ! local version of kd (normal counter)
        character(len=max_str_ln) :: plot_title                                 ! title for plots
        real(dp), allocatable :: ranges(:,:)                                    ! the range of the plots
        real(dp), allocatable :: theta_loc(:)                                   ! local version of theta
        real(dp) :: q_saf_loc                                                   ! safety factor at current normal point
        real(dp) :: norm_pt                                                     ! normal point
        
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
            
            ! set up n_norm and inc_norm
            n_norm = min(n_r_X,max_n_norm)
            inc_norm = max(1.0_dp,(n_r_X-1.0_dp)/(max_n_norm-1))
            
            ! set up  theta, depending on whether it has  been provided as input
            ! or not
            if (present(theta)) then
                allocate(theta_loc(n_theta))
                theta_loc = [(theta(1) + (jd-1.0_dp)/(n_theta-1.0_dp) * &
                    &(theta(2)-theta(1)),jd=1,n_theta)]
            else
                n_theta = 1
                allocate(theta_loc(n_theta))
                theta_loc = 0.0_dp
            end if
            
            ! set up x-axis values
            allocate(x_plot(n_norm,n_theta,n_m_X,n_n_t*n_t))
            do id = 1,n_norm
                x_plot(id,:,:,:) = min_r_X + (id-1.0)/(n_norm-1)*(max_r_X-min_r_X)
            end do
            
            ! set up y-axis values
            allocate(y_plot(n_norm,n_theta,n_m_X,n_n_t*n_t))
            if (present(theta)) then
                do id = 1,n_theta
                    y_plot(:,id,:,:) = theta_loc(id)
                end do
            else
                y_plot = 0.0_dp
            end if
            
            ! set up z-axis values (not same for all time steps)
            ! allocate variables
            allocate(z_plot(n_norm,n_theta,n_m_X,n_n_t*n_t))
            allocate(max_of_modes(n_m_X))
            allocate(max_of_modes_x(n_m_X))
            allocate(max_of_modes_y(n_m_X))
            max_of_modes = 0.0_dp
            max_of_modes_x = 0.0_dp
            max_of_modes_y = 0.0_dp
            
            z_plot = 10
            ! plot for all time steps the Eigenfunction at theta = 0
            do id = 1,n_n_t*n_t
                ! loop over all normal points of all modes
                kd_loc = 1.0_dp
                do kd = 1,n_norm
                    ! loop over all modes
                    do ld = 1,n_m_X
                        do jd = 1,n_theta
                            call dis2con(kd,[1,n_norm],norm_pt,[min_r_X,max_r_X])
                            ierr = calc_interp(q_saf_full(:,0),[1,n_r_eq],&
                                &q_saf_loc,norm_pt)
                            CHCKERR('')
                            z_plot(kd,jd,ld,id) = realpart(&
                                &ser_X_vec(ld,nint(kd_loc)) * &
                                &exp(omega/abs(omega)*(id-1.0_dp)/(n_t*n_n_t-1)&
                                & + iu*(n_X*q_saf_loc-m_X(ld))*theta_loc(jd)))
                            if (id.eq.1 .or. id.eq.n_t+1) then
                                if (z_plot(kd,jd,ld,id).gt.max_of_modes(ld)) then
                                    max_of_modes(ld) = z_plot(kd,jd,ld,id)
                                    max_of_modes_x(ld) = x_plot(kd,jd,ld,id)
                                    max_of_modes_y(ld) = y_plot(kd,jd,ld,id)
                                end if
                            end if
                        end do
                    end do
                    kd_loc = kd_loc + inc_norm
                end do
                if (id.eq.1 .or. id.eq.n_t+1) then
                    !! plot places of maximum
                    if (id.eq.1) then
                        !call writo('Information about the maximum of the &
                            !&modes for time t=0')
                        plot_title = 'JOB '//trim(i2str(alpha_job_nr))//&
                            &' - Eigenvalue '//trim(i2str(X_id))//&
                            &' - Eigenvector (t = 0)'
                    else if (id.eq.n_t+1) then
                        !call writo('Information about the maximum of the &
                            !&modes for time t=T_0/4')
                        plot_title = 'JOB '//trim(i2str(alpha_job_nr))//&
                            &' - Eigenvalue '//trim(i2str(X_id))//&
                            &' - Eigenvector (t = 0.25 T_0)'
                    end if
                    !call print_GP_2D('maximum of the modes','',max_of_modes)
                    !call print_GP_2D('normal flux surface where max. of modes &
                        !&occurs','',max_of_modes_x)
                    !call print_GP_2D('parallel point where max. of modes &
                        !&occurs','',max_of_modes_y)
                    ! plot the modes on screen
                    if (present(theta)) then
                        call print_GP_3D(plot_title,'Eigenvector_'//&
                            &trim(i2str(alpha_job_nr))//'.dat',z_plot(:,:,:,id)&
                            &,y=y_plot(:,:,:,id),x=x_plot(:,:,:,id),draw=.true.)
                    else
                        call print_GP_2D(plot_title,'Eigenvector_'//&
                            &trim(i2str(alpha_job_nr))//'.dat',z_plot(:,1,:,id)&
                            &,x=x_plot(:,1,:,id),draw=.true.)
                    end if
                    ! plot them in file as well
                    call draw_GP(plot_title,'Eigenvector_'//&
                        &trim(i2str(alpha_job_nr))//'.dat',n_m_X,&
                        &.not.present(theta),.false.)
                end if
            end do
            
            deallocate(max_of_modes)
            deallocate(max_of_modes_x)
            
            ! allocate ranges
            if (present(theta)) then
                allocate(ranges(3,2))
            else
                allocate(ranges(2,2))
            end if
            
            if (n_m_X.gt.1) then
                ! plot the individual modes in time
                do ld = 1,min(n_m_X,max_n_plots)
                    ! set plot_title
                    plot_title = 'JOB '//trim(i2str(alpha_job_nr))//&
                        &' - mode '//trim(i2str(ld))//' of '//trim(i2str(n_m_X))
                    if (present(theta)) then
                        ranges(1,:) = [minval(x_plot(:,:,ld,:)),&
                            &maxval(x_plot(:,:,ld,:))]
                        ranges(2,:) = [minval(y_plot(:,:,ld,:)),&
                            &maxval(y_plot(:,:,ld,:))]
                        ranges(3,:) = [minval(z_plot(:,:,ld,:)),&
                            &maxval(z_plot(:,:,ld,:))]
                        call print_GP_3D(trim(plot_title),'temp_mode_'//&
                            &trim(i2str(alpha_job_nr)),z_plot(:,:,ld,:),&
                            &y=y_plot(:,:,ld,:),x=x_plot(:,:,ld,:),draw=.false.)
                    else
                        ranges(1,:) = [minval(x_plot(:,:,ld,:)),&
                            &maxval(x_plot(:,:,ld,:))]
                        ranges(2,:) = [minval(z_plot(:,:,ld,:)),&
                            &maxval(z_plot(:,:,ld,:))]
                        call print_GP_2D(trim(plot_title),'temp_mode_'//&
                            &trim(i2str(alpha_job_nr)),z_plot(:,1,ld,:),&
                            &x=x_plot(:,1,ld,:),draw=.false.)
                    end if
                    call draw_GP_animated(trim(plot_title),'temp_mode_'//&
                        &trim(i2str(alpha_job_nr)),n_n_t*n_t,.false.,ranges)
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
            do ld = 2,n_m_X
                z_plot(:,:,1,:) = z_plot(:,:,1,:) + z_plot(:,:,ld,:)
            end do
            plot_title = 'JOB '//trim(i2str(alpha_job_nr))//' - global mode'
            if (present(theta)) then
                ranges(1,:) = [minval(x_plot(:,:,1,:)),maxval(x_plot(:,:,1,:))]
                ranges(2,:) = [minval(y_plot(:,:,1,:)),maxval(y_plot(:,:,1,:))]
                ranges(3,:) = [minval(z_plot(:,:,1,:)),maxval(z_plot(:,:,1,:))]
                call print_GP_3D(trim(plot_title),'temp_mode_'//&
                    &trim(i2str(alpha_job_nr)),z_plot(:,:,1,:),&
                    &y=y_plot(:,:,1,:),x=x_plot(:,:,1,:),draw=.false.)
            else
                ranges(1,:) = [minval(x_plot(:,:,1,:)),maxval(x_plot(:,:,1,:))]
                ranges(2,:) = [minval(z_plot(:,:,1,:)),maxval(z_plot(:,:,1,:))]
                call print_GP_2D(trim(plot_title),'temp_mode_'//&
                    &trim(i2str(alpha_job_nr)),z_plot(:,1,1,:),&
                    &x=x_plot(:,1,1,:),draw=.false.)
            end if
            call draw_GP_animated(trim(plot_title),'temp_mode_'//&
                &trim(i2str(alpha_job_nr)),n_n_t*n_t,.not.present(theta),ranges)
            
            ! deallocate variables
            deallocate(ranges)
            deallocate(ser_X_vec)
            deallocate(x_plot)
            deallocate(y_plot)
            deallocate(z_plot)
        end if
        
        call lvl_ud(-1)
    end function plot_X_vec
end module
