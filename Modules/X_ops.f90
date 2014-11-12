!------------------------------------------------------------------------------!
!   Calculate the matrix  elements due to the Plasma of  the system of equations
!   that has to be solved as described in [ADD REF]
!------------------------------------------------------------------------------!
module X_ops
#include <PB3D_macros.h>
    use num_vars, only: dp, iu, max_str_ln, pi
    use message_ops, only: lvl_ud, writo, print_ar_2
    use output_ops, only: print_GP_2D, print_GP_3D
    use str_ops, only: i2str, r2strt, r2str

    implicit none
    private
    public prepare_X, solve_EV_system, plot_X_vec

contains
    ! prepare the matrix elements by calculating KV and PV, which then will have
    ! to be integrated, with a complex exponential weighting function
    subroutine prepare_X
        use X_vars, only: init_X, calc_PV, calc_KV, calc_U, calc_extra, &
            &calc_V_int, dealloc_X, &
            &PV0, PV1, PV2, KV0, KV1, KV2, PV_int, KV_int
        use num_vars, only: use_pol_flux
        
        ! local variables
        character(len=5) :: ang_par_F_name                                      ! error message
        
        ! set up ang_par_F_name
        if (use_pol_flux) then
            ang_par_F_name = 'theta'
        else
            ang_par_F_name = 'zeta'
        end if
        
        call writo('Calculating table of field line averages &
            &<V e^i[(k-m)'//trim(ang_par_F_name)//']>')
        
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
        
        ! Calculate PV_int = <PVi e^(k-m)ang_par_F>
        call writo('Taking field average of PV')
        call lvl_ud(1)
        call calc_V_int(PV0,PV_int(:,:,:,1))
        call calc_V_int(PV1,PV_int(:,:,:,2))
        call calc_V_int(PV2,PV_int(:,:,:,3))
        call lvl_ud(-1)
        
        ! Calculate KV_int = <KVi e^(k-m)ang_par_F>
        call writo('Taking field average of KV')
        call lvl_ud(1)
        call calc_V_int(KV0,KV_int(:,:,:,1))
        call calc_V_int(KV1,KV_int(:,:,:,2))
        call calc_V_int(KV2,KV_int(:,:,:,3))
        call lvl_ud(-1)
        
        ! deallocate equilibrium variables
        call writo('deallocating unused variables')
        call dealloc_X
        
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
    
    ! Plots an Eigenfunction
    ! [MPI] Collective call
    integer function plot_X_vec(X_vec,X_val,X_id,job_id,par_range_F) &
        &result(ierr)
        use num_vars, only: output_style
        
        character(*), parameter :: rout_name = 'plot_X_vec'
        
        ! input / output
        complex(dp), intent(in) :: X_vec(:,:)                                   ! MPI Eigenvector
        complex(dp), intent(in) :: X_val                                        ! Eigenvalue
        integer, intent(in) :: X_id                                             ! nr. of Eigenvalue (for output name)
        integer, intent(in) :: job_id                                           ! nr. of alpha job
        real(dp), intent(in), optional :: par_range_F(2)                        ! parallel range [2pi]
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! local variables
        integer :: n_t(2)                                                       ! nr. of time steps in quarter period, nr. of quarter periods
        complex(dp) :: omega                                                    ! sqrt of Eigenvalue
        
        ! initialize ierr
        ierr = 0
        
        ! if omega is  predominantly real, the Eigenfunction  explodes, so limit
        ! n_t(2)  to 1.  If  it is  predominantly  imaginary, the  Eigenfunction
        ! oscillates, so choose n_t(2) = 8 for 2 whole periods
        omega = sqrt(X_val)
        if (abs(realpart(omega)/imagpart(omega)).gt.1) then
            n_t(2) = 1                                                          ! 1 quarter period
            n_t(1) = 10                                                         ! 10 points per quarter period
        else
            n_t(2) = 8                                                          ! 8 quarter periods
            n_t(1) = 5                                                          ! 5 points per quarter period
        end if
        
        ! delegate to child routines
        select case(output_style)
            case(1)                                                             ! GNUPlot output
                ierr = plot_X_vec_GP(X_vec,omega,X_id,job_id,n_t,par_range_F)
                CHCKERR('')
            case(2)                                                             ! HDF5 output
                ierr = plot_X_vec_HDF5(X_vec,X_val,X_id,job_id,n_t,.false.)
                CHCKERR('')
            case default
                err_msg = 'No style associated with '//trim(i2str(output_style))
                ierr = 1
                CHCKERR(err_msg)
        end select
    end function plot_X_vec
    
    ! GNUPlot version of plot_X_vec
    ! plots an  eigenfunction in either 1D  (normal coord. in perturb.  grid) or
    ! 2D  (normal  coord. &  par  coord.),  depending  on whether  the  variable
    ! par_range_F is provided The individual modes are  plotted at t = 0 and t =
    ! (pi/2)/omega and  an animated plot  is produced as  well in a  .gif output
    ! file
    integer function plot_X_vec_GP(X_vec,omega,X_id,job_id,n_t,par_range_F) &
        &result(ierr)
        use num_vars, only: grp_rank, min_r_X, max_r_X, max_n_plots, &
            &grp_n_procs, use_pol_flux
        use output_ops, only: draw_GP, draw_GP_animated, merge_GP
        use X_vars, only: n_r_X, size_X, n_X, m_X, grp_r_X, grp_min_r_X
        use VMEC_vars, only: n_r_eq, VMEC_use_pol_flux
        use MPI_ops, only: get_ser_X_vec, get_ghost_X_vec, wait_MPI
        use eq_vars, only: q_saf_FD, rot_t_FD, grp_n_r_eq, grp_min_r_eq, &
            &max_flux, max_flux_VMEC, flux_p_FD, flux_t_FD
        use utilities, only: dis2con, con2dis, interp_fun_1D
        
        character(*), parameter :: rout_name = 'plot_X_vec_GP'
        
        ! input / output
        complex(dp), intent(in) :: X_vec(:,:)                                   ! MPI Eigenvector
        complex(dp), intent(in) :: omega                                        ! sqrt of Eigenvalue
        integer, intent(in) :: X_id                                             ! nr. of Eigenvalue (for output name)
        integer, intent(in) :: job_id                                           ! nr. of alpha job
        integer, intent(in) :: n_t(2)                                           ! nr. of time steps in quarter period, nr. of quarter periods
        real(dp), intent(in), optional :: par_range_F(2)                        ! parallel range [2pi]
        
        ! local variables
        complex(dp), allocatable :: X_vec_extended(:,:)                         ! MPI Eigenvector extended with assymetric ghost region
        real(dp), allocatable :: x_plot(:,:,:,:)                                ! x values of plot, parallel version
        real(dp), allocatable :: y_plot(:,:,:,:)                                ! y values of plot, parallel version
        real(dp), allocatable :: z_plot(:,:,:,:)                                ! z values of plot, parallel version
        real(dp), allocatable :: z_magn_plot(:,:,:)                             ! z values of plot of magnitude, parallel version
        integer :: id, jd, kd, ld                                               ! counter
        integer :: n_par_F = 50                                                 ! how many parallel points
        real(dp) :: kd_loc_eq                                                   ! kd_loc in equilibrium grid
        real(dp) :: kd_loc_i                                                    ! unrounded index corresponding to kd_loc in tables
        character(len=max_str_ln) :: plot_title                                 ! title for plots
        character(len=max_str_ln) :: file_name                                  ! file name for plots
        real(dp), allocatable :: ang_par_F_loc(:)                               ! local version of ang_par_F
        complex(dp), allocatable :: X_vec_interp(:)                             ! X_vec at an interpolated normal position
        complex(dp), allocatable :: first_X_vec_interp(:)                       ! first X_vec_interp
        real(dp), allocatable :: fac_n(:), fac_m(:)                             ! multiplicative factors for n and m
        real(dp) :: fac_n_interp, fac_m_interp                                  ! fac_n and fac_m at interpolated normal position
        real(dp) :: first_fac_n_interp, first_fac_m_interp                      ! first fac_n_interp and fac_m_interp
        character(len=max_str_ln), allocatable :: file_names(:)                 ! names of file of plots of different ranks
        real(dp), pointer :: flux(:), flux_VMEC(:)                              ! either pol. or tor. flux
        real(dp), allocatable :: r_plot(:)                                      ! local normal values at which to interpolate for the plot
        integer :: n_norm                                                       ! max nr of normal points to plot
        integer :: grp_n_norm                                                   ! nr. of normal points in plot for this rank
        complex(dp) :: fac_i_interp(1)                                          ! complex copy of fac_n_interp or fac_m_interp
        
        ! initialize ierr
        ierr = 0
        
        call lvl_ud(1)
        
        ! 1. initialize some quantities
        
        ! set up fac_n and fac_m
        allocate(fac_n(grp_n_r_eq),fac_m(grp_n_r_eq))
        if (use_pol_flux) then
            fac_n = q_saf_FD(:,0)
            fac_m = 1.0_dp
        else
            fac_n = 1.0_dp
            fac_m = rot_t_FD(:,0)
        end if
        
        ! set up the angular variable ang_par_F_loc, depending on whether it has
        ! been provided as input or not
        if (present(par_range_F)) then
            allocate(ang_par_F_loc(n_par_F))
            ang_par_F_loc = &
                &[(par_range_F(1) + (jd-1.0_dp)/(n_par_F-1.0_dp) * &
                &(par_range_F(2)-par_range_F(1)),jd=1,n_par_F)]
        else
            n_par_F = 1
            allocate(ang_par_F_loc(n_par_F))
            ang_par_F_loc = 0.0_dp
        end if
        
        ! 2. calculate the x, y and z parts of the plot in parallel
        
        ! calculate starting and ending normal index of array
        call calc_r_plot(r_plot,n_norm,grp_n_norm,.true.)
        
        ! allocate variables
        if (grp_rank.lt.grp_n_procs-1) then
            allocate(X_vec_extended(size(X_vec,1),size(X_vec,2)+1))
        else
            allocate(X_vec_extended(size(X_vec,1),size(X_vec,2)))
        end if
        allocate(x_plot(grp_n_norm,n_par_F,size_X,product(n_t)))
        allocate(y_plot(grp_n_norm,n_par_F,size_X,product(n_t)))
        allocate(z_plot(grp_n_norm,n_par_F,size_X,product(n_t)))
        allocate(z_magn_plot(grp_n_norm,n_par_F,product(n_t)))
        allocate(X_vec_interp(size_X))
        allocate(first_X_vec_interp(size_X))
        
        ! set up extended X_vec
        X_vec_extended(:,1:size(X_vec,2)) = X_vec
        ierr = get_ghost_X_vec(X_vec_extended,1)
        CHCKERR('')
        
        ! set up y-axis in parallel
        if (present(par_range_F)) then
            do id = 1,n_par_F
                y_plot(:,id,:,:) = ang_par_F_loc(id)
            end do
        else
            y_plot = 0.0_dp
        end if
        
        ! set up flux and flux_VMEC
        if (use_pol_flux) then
            flux => flux_p_FD(:,0)
        else
            flux => flux_t_FD(:,0)
        end if
        if (VMEC_use_pol_flux) then
            flux_VMEC => flux_p_FD(:,0)
        else
            flux_VMEC => flux_t_FD(:,0)
        end if
        
        ! set up x-axis and z-axis values in parallel
        z_magn_plot = 0.0_dp
        ! loop over all time steps
        do id = 1,product(n_t)
            ! loop over all normal points
            do kd = 1,grp_n_norm
                ! set up x_plot
                x_plot(kd,:,:,:) = r_plot(kd)
                
                ! set up the interpolated values of X_vec, fac_n and fac_m
                if (kd.lt.grp_n_norm .and. grp_rank+1.lt.grp_n_procs &
                    &.or. grp_rank+1.eq.grp_n_procs) then                       ! last point is treated differently
                    ! set up interp. X_vec (tabulated in perturbation grid)
                    call con2dis(r_plot(kd),[min_r_X,max_r_X],kd_loc_i,&
                        &[1,n_r_X])
                    kd_loc_i = kd_loc_i - grp_min_r_X + 1                       ! include offset of perturbation tables
                    X_vec_interp = X_vec_extended(:,floor(kd_loc_i)) + &
                        &(kd_loc_i-floor(kd_loc_i)) * &
                        &(X_vec_extended(:,ceiling(kd_loc_i))-&
                        &X_vec_extended(:,floor(kd_loc_i)))
                    
                    ! set up  interp. fac_n and fac_m  (tabulated in equilibrium
                    ! grid)
                    ierr = interp_fun_1D(kd_loc_eq,flux_VMEC/max_flux_VMEC,&
                        &r_plot(kd),flux/max_flux)
                    CHCKERR('')
                    call con2dis(kd_loc_eq,[0._dp,1._dp],kd_loc_i,[1,n_r_eq])
                    kd_loc_i = kd_loc_i - grp_min_r_eq + 1                      ! include offset of equilibrium tables
                    fac_n_interp = fac_n(floor(kd_loc_i)) + &
                        &(kd_loc_i-floor(kd_loc_i)) * &
                        &(fac_n(ceiling(kd_loc_i))-fac_n(floor(kd_loc_i)))
                    fac_m_interp = fac_m(floor(kd_loc_i)) + &
                        &(kd_loc_i-floor(kd_loc_i)) * &
                        &(fac_m(ceiling(kd_loc_i))-fac_m(floor(kd_loc_i)))
                end if
                
                ! save first interpolated X_vec, fac_n and fac_m
                if (kd.eq.1) then
                    first_X_vec_interp = X_vec_interp
                    first_fac_n_interp = fac_n_interp
                    first_fac_m_interp = fac_m_interp
                end if
                
                ! exchange X_vec_interp, fac_n and fac_m at last points
                if (kd.eq.grp_n_norm) then
                    ! X_vec_interp
                    call exch_ghost_values(first_X_vec_interp)
                    if (grp_rank+1.lt.grp_n_procs) &
                        &X_vec_interp = first_X_vec_interp
                    
                    ! fac_n_interp
                    fac_i_interp = first_fac_n_interp
                    call exch_ghost_values(fac_i_interp)
                    if (grp_rank+1.lt.grp_n_procs) &
                        &fac_n_interp = realpart(fac_i_interp(1))
                    
                    ! fac_m_interp
                    fac_i_interp = first_fac_m_interp
                    call exch_ghost_values(fac_i_interp)
                    if (grp_rank+1.lt.grp_n_procs) &
                        &fac_m_interp = realpart(fac_i_interp(1))
                end if
                
                ! loop over all modes and all parallel points to set up z_plot
                do ld = 1,size_X
                    do jd = 1,n_par_F
                        z_plot(kd,jd,ld,id) = &
                            &realpart(X_vec_interp(ld) * &
                            &exp(omega/abs(omega)*0.5*pi*(id-1.0_dp)/(n_t(1)-1)&
                            &    + iu*(n_X(ld)*fac_n_interp-&
                            &          m_X(ld)*fac_m_interp)&
                            &    * ang_par_F_loc(jd)))
                    end do
                    z_magn_plot(kd,:,id) = z_magn_plot(kd,:,id) + &
                        &z_plot(kd,:,ld,id)
                end do
            end do
        end do
        
        ! 3. for  every mode, write plot  output data in data files  and plot by
        ! group master
        if (size_X.gt.1) then
            ! plot the individual modes in time
            do ld = 1,min(size_X,max_n_plots)
                ! set plot_title and file name and write data
                plot_title = 'job '//trim(i2str(job_id))//' - EV '//&
                    &trim(i2str(X_id))//' - mode '//trim(i2str(ld))//' of '//&
                    &trim(i2str(size_X))
                file_name = 'mode_'//trim(i2str(X_id))//'_'//&
                    &trim(i2str(ld))//'.dat'
                if (present(par_range_F)) then
                    call print_GP_3D(trim(plot_title),trim(file_name)//'_'//&
                        &trim(i2str(grp_rank)),z_plot(:,:,ld,:),&
                        &y=y_plot(:,:,ld,:),x=x_plot(:,:,ld,:),draw=.false.)
                else
                    call print_GP_2D(trim(plot_title),trim(file_name)//'_'//&
                        &trim(i2str(grp_rank)),z_plot(:,1,ld,:),&
                        &x=x_plot(:,1,ld,:),draw=.false.)
                end if
                    
                ! wait for all processes
                ierr = wait_MPI()
                CHCKERR('')
                
                ! plot by group master
                if (grp_rank.eq.0) then
                    ! set up file names in array
                    allocate(file_names(grp_n_procs))
                    do id = 1,grp_n_procs
                        file_names(id) = trim(file_name)//'_'//trim(i2str(id-1))
                    end do
                    
                    ! merge files
                    call merge_GP(file_names,file_name,delete=.true.)
                    
                    ! draw animation
                    call draw_GP_animated(trim(plot_title),file_name,&
                        &product(n_t),.false.)
                    
                    ! deallocate
                    deallocate(file_names)
                end if
            end do
            if (size_X.eq.max_n_plots+1) then
                call writo('Skipping plot for mode '//&
                    &trim(i2str(max_n_plots+1)))
            else if (size_X.gt.max_n_plots+1) then
                call writo('Skipping plots for modes '//&
                    &trim(i2str(max_n_plots+1))//'..'//trim(i2str(size_X)))
            end if
        end if
        
        ! 4. do the same for the global mode
        ! set plot_title and file name and write data
        plot_title = 'job '//trim(i2str(job_id))//' - EV '//&
            &trim(i2str(X_id))//' - global mode'
        file_name = 'global_mode_'//trim(i2str(X_id))//'.dat'
        if (present(par_range_F)) then
            call print_GP_3D(trim(plot_title),trim(file_name)//&
                &'_'//trim(i2str(grp_rank)),&
                &z_magn_plot(:,:,:),y=y_plot(:,:,1,:),&
                &x=x_plot(:,:,1,:),draw=.false.)
        else
            call print_GP_2D(trim(plot_title),trim(file_name)//&
                &'_'//trim(i2str(grp_rank)),&
                &z_magn_plot(:,1,:),x=x_plot(:,1,1,:),&
                &draw=.false.)
        end if
        
        ! plot by group master
        if (grp_rank.eq.0) then
            ! set up file names in file
            allocate(file_names(grp_n_procs))
            do id = 1,grp_n_procs
                file_names(id) = trim(file_name)//'_'//trim(i2str(id-1))
            end do
            
            ! merge files
            call merge_GP(file_names,file_name,delete=.true.)
            
            ! draw animation
            call draw_GP_animated(trim(plot_title),file_name,product(n_t),&
                &.false.)
            
            ! deallocate
            deallocate(file_names)
        end if
        
        call lvl_ud(-1)
    contains
        ! calculates the array  of normal points in the  perturbation grid where
        ! to interpolate  to produce the plot.  This assumes that the  vector to
        ! be  plot  has  an  assymetric  ghost  region  of  size  one,  so  that
        ! the  interpolated  plot points  can  be  calculated in  every  process
        ! independently
        ! Additionally,  there can be  an extra  assymetric ghost region  in the
        ! plot grid so that there is overlap (indicated by extra_ghost)
        subroutine calc_r_plot(r_plot,n_norm,grp_n_norm,extra_ghost)
            ! input / output
            real(dp), intent(inout), allocatable :: r_plot(:)                   ! local normal values at which to interpolate for the plot
            integer, intent(inout) :: grp_n_norm                                ! nr. of normal points for this rank
            integer, intent(inout) :: n_norm                                    ! total nr. of normal points
            logical, intent(in), optional :: extra_ghost                        ! .true. if extra ghost region requested (e.g. for GNUPlot)
            
            ! local variables
            integer :: start_i_norm, end_i_norm                                 ! starting and ending index of normal points for this rank
            integer :: max_n_norm = 50                                          ! nr. of normal points
            real(dp) :: inc_norm                                                ! normal increment for plots
            integer :: kd                                                       ! counter
            
            ! set  up how  many  normal  points for  this  rank  n_norm and  the
            ! increment between normal points inc_norm
            ! (taking all the normal points could be way too many)
            n_norm = min(n_r_X,max_n_norm)
            inc_norm = (max_r_X-min_r_X)/(n_norm-1)
            
            ! calculate the  starting and ending index,  first determining which
            ! is the range where normal interpolation can take place
            ! from:
            !   gpr_r_X(1) .le. (min_r_X+(k-1)*inc_norm) .le. grp_r_X(grp_n_r_X)
            ! with k an integer
            start_i_norm = ceiling((grp_r_X(1)-min_r_X)/inc_norm + 1)
            end_i_norm = floor((grp_r_X(size(grp_r_X))-min_r_X)/inc_norm + 1)   ! size(grp_r_X) can be grp_n_r_X or grp_n_r_X + 1 (ghost region in X grid)
            
            ! increment end_i_norm with one if not  last rank, so that there can
            ! be overlap in the plots
            if (present(extra_ghost)) then
                if (extra_ghost) then
                    if (grp_rank.lt.grp_n_procs-1) end_i_norm = end_i_norm + 1
                end if
            end if
            
            ! set grp_n_norm
            grp_n_norm = end_i_norm-start_i_norm+1
            
            ! allocate r_plot with the correct size
            allocate(r_plot(grp_n_norm)); 
            
            ! set r_plot
            r_plot = [(min_r_X + (start_i_norm+kd-2)*inc_norm,kd=1,grp_n_norm)]
            if (r_plot(grp_n_norm).gt.1.0_dp) r_plot(grp_n_norm) = 1.0_dp
        end subroutine calc_r_plot
        
        ! exhange the ghost values of a vector
        subroutine exch_ghost_values(vec)
            ! input / output
            complex(dp), intent(inout) :: vec(:)                                ! input vector, later overwritten with exchanged value
            
            ! local variables
            complex(dp), allocatable :: exch_vec(:,:)                           ! to exchange the input vector
            integer :: n_modes                                                  ! number of modes to exchange
            
            ! set up n_modes
            n_modes = size(vec)
            
            ! allocate the exchange vectors
            if (grp_rank.eq.0) then
                allocate(exch_vec(n_modes,1))
            else if (grp_rank+1.eq.grp_n_procs) then
                allocate(exch_vec(n_modes,1))
                exch_vec(:,1) = vec                                             ! fill with vec
            else
                allocate(exch_vec(n_modes,2))
                exch_vec(:,1) = vec                                             ! fill with vec
            end if
            
            ! exchange the last interpolated X_vec
            ierr = get_ghost_X_vec(exch_vec,1)
            CHCKERR('')
            
            ! copy value to new interpolated X_vec
            vec = exch_vec(:,size(exch_vec,2))
        end subroutine exch_ghost_values
    end function plot_X_vec_GP
    
    ! HDF5 version of plot_X_vec
    ! Plots an  eigenfunction in  3D real  space by determining  a grid  in VMEC
    ! coordinates that covers  the entire device (But they are  tabulated in the
    ! perturbation normal  grid) and then  calculating the real  perturbation on
    ! that grid  by inverting the Fourier  transform that defined the  modes for
    ! which are  solved.
    ! Alternatively, by using the optional argument "follow_B", the modes can be
    ! plot along  the magnetic  field lines  which have been  used to  solve the
    ! system of equations. (TO BE IMPLEMENTED)
    integer function plot_X_vec_HDF5(X_vec,omega,X_id,job_id,n_t,follow_B) &
        &result(ierr)
        use X_vars, only: grp_min_r_X, grp_n_r_X, n_r_X, grp_r_X, grp_n_r_X, &
            &size_X, n_X, m_X
        use output_ops, only: print_HDF5_3D
        use eq_vars, only: calc_XYZ_grid
        use HDF5_vars, only: open_HDF5_file, print_HDF5_top, print_HDF5_geom, &
            &print_HDF5_3D_data_item, print_HDF5_grid, add_HDF5_item, &
            &close_HDF5_file, print_HDF5_att, &
            &HDF5_file_type, XML_str_type
        
        character(*), parameter :: rout_name = 'plot_X_vec_HDF5'
        
        ! input / output
        complex(dp), intent(in) :: X_vec(:,:)                                   ! MPI Eigenvector
        complex(dp), intent(in) :: omega                                        ! sqrt of Eigenvalue
        integer, intent(in) :: X_id                                             ! nr. of Eigenvalue (for output name)
        integer, intent(in) :: job_id                                           ! nr. of alpha job
        integer, intent(in) :: n_t(2)                                           ! nr. of time steps in quarter period, nr. of quarter periods
        logical, intent(in), optional :: follow_B                               ! .true. if to be plot only along B
        
        ! local variables
        integer :: id, kd, ld                                                   ! counters
        integer :: n_theta_plot = 80                                            ! nr. of poloidal points in plot
        integer :: n_zeta_plot = 160                                            ! nr. of toroidal points in plot
        real(dp), allocatable :: theta_plot(:,:,:), zeta_plot(:,:,:)            ! theta_V and zeta_V for flux surface plot
        real(dp), allocatable :: x_plot(:,:,:)                                  ! x values of plot
        real(dp), allocatable :: y_plot(:,:,:)                                  ! y values of plot
        real(dp), allocatable :: z_plot(:,:,:)                                  ! z values of plot
        real(dp), allocatable :: l_plot(:,:,:)                                  ! lambda values of plot
        real(dp), allocatable :: f_plot(:,:,:)                                  ! the function to plot
        logical :: follow_B_loc                                                 ! local copy of follow_B
        integer :: tot_dim(3), grp_dim(3), grp_offset(3)                        ! total dimensions, group dimensions and group offset
        type(HDF5_file_type) :: file_info                                       ! HDF5 file info
        type(XML_str_type) :: time_col_grid                                     ! grid with time collection
        type(XML_str_type), allocatable :: grids(:)                             ! the grids in the time collection
        type(XML_str_type) :: top                                               ! topology
        type(XML_str_type) :: XYZ(3)                                            ! data items for geometry
        type(XML_str_type) :: geom                                              ! geometry
        type(XML_str_type) :: att(1)                                            ! attribute
        character(len=max_str_ln) :: var_name                                   ! name of variable that is plot
        
        ! initialize ierr
        ierr = 0
        
        call lvl_ud(1)
        
        ! user output
        call writo('Starting the plot')
        
        ! set up local follow_B
        follow_B_loc = .false.
        if (present(follow_B)) then
            if (follow_B) follow_B_loc = follow_B
        end if
        
        ! set up total dimensions, group dimensions and group offset
        tot_dim = [n_theta_plot,n_zeta_plot,n_r_X]
        grp_dim = [n_theta_plot,n_zeta_plot,grp_n_r_X]
        grp_offset = [0,0,grp_min_r_X-1]
        
        ! set up var_name
        var_name = 'Solution vector X_vec'
        
        ! initialize theta_plot and zeta_plot
        if (follow_B) then
            ierr = 1
            CHCKERR('NOT YET IMPLEMENTED')
            !!!!! SHOULD BE EQUAL TO THE GRID PLOT, WITH ADDITIONALLY THE SOLUTION VECTOR !!!!
        else
            allocate(theta_plot(n_theta_plot,n_zeta_plot,grp_n_r_X))
            do id = 1,n_theta_plot
                theta_plot(id,:,:) = (id-1.0_dp)*2*pi/(n_theta_plot-1.0_dp)
            end do
            allocate(zeta_plot(n_theta_plot,n_zeta_plot,grp_n_r_X))
            do id = 1,n_zeta_plot
                zeta_plot(:,id,:) = (id-1.0_dp)*2*pi/(n_zeta_plot-1.0_dp)
            end do
        end if
        
        ! calculate X,Y and Z using the VMEC theta_plot and zeta_plot, tabulated
        ! in the perturbation normal grid indicated by grp_r_X
        ierr = calc_XYZ_grid(theta_plot,zeta_plot,&
            &[grp_r_X(1),grp_r_X(grp_n_r_X)],.false.,&
            &x_plot,y_plot,z_plot,l_plot)
        CHCKERR('')
        
        !call print_GP_3D('test_plot','',z_plot,X=x_plot,Y=y_plot)               ! for testing
        
        ! calculate the function to plot: global mode
        allocate(f_plot(n_theta_plot,n_zeta_plot,grp_n_r_X))
        
        ! open HDF5 file
        ierr = open_HDF5_file(file_info,'X_vec_'//trim(i2str(job_id)),&
            &'Job '//trim(i2str(job_id))//' - Solution vector X_vec &
            &for Eigenvalue '//trim(i2str(X_id)))
        CHCKERR('')
        
        ! print topology
        call print_HDF5_top(top,2,tot_dim)
        
        ! add topology to HDF5 file and reset it
        call add_HDF5_item(file_info,top,reset=.true.)
        
        ! print data item for X
        ierr = print_HDF5_3D_data_item(XYZ(1),file_info,&
            &'X_'//trim(i2str(id)),x_plot,tot_dim,grp_dim=grp_dim,&
            &grp_offset=grp_offset)
        CHCKERR('')
        
        ! print data item for Y
        ierr = print_HDF5_3D_data_item(XYZ(2),file_info,&
            &'Y_'//trim(i2str(id)),y_plot,tot_dim,grp_dim=grp_dim,&
            &grp_offset=grp_offset)
        CHCKERR('')
        
        ! print data item for Z
        ierr = print_HDF5_3D_data_item(XYZ(3),file_info,&
            &'Z_'//trim(i2str(id)),z_plot,tot_dim,grp_dim=grp_dim,&
            &grp_offset=grp_offset)
        CHCKERR('')
        
        ! print geometry with X, Y and Z data item
        call print_HDF5_geom(geom,2,XYZ,.true.)
        
        ! add geometry to HDF5 file and reset it
        call add_HDF5_item(file_info,geom,reset=.true.)
        
        ! create grid for time collection
        allocate(grids(product(n_t)))
        
        ! user output
        call writo('Calculating plot for '//trim(i2str(product(n_t)))//&
            &' time steps')
        
        call lvl_ud(1)
            
        ! For  each time step, produce the 3D plot
        do id = 1,product(n_t)
            ! print data item for plot variable f_plot
            f_plot = 0.0_dp
            ! for all normal points
            do kd = 1,grp_n_r_X
                ! for all modes
                do ld = 1,size_X
                    ! Need to  translate from VMEC  coordinates (theta_V,zeta_V)
                    ! to  flux  coordinates  (theta_F,zeta_F)  for  the  inverse
                    ! Fourier transform of the  Eigenvalue Fourier modes to real
                    ! space
                    f_plot(:,:,kd) = f_plot(:,:,kd) + realpart(X_vec(ld,kd) * &
                        &exp(omega/abs(omega)*0.5*pi*(id-1.0_dp)/(n_t(1)-1) + &
                        &iu * (n_X(ld)*(-zeta_plot(:,:,kd)) - &
                        &m_X(ld)*(theta_plot(:,:,kd)+l_plot(:,:,kd)))))
                end do
            end do
            
            ierr = print_HDF5_3D_data_item(XYZ(1),file_info,'var_'//&
                &trim(i2str(id)),f_plot,tot_dim,grp_dim=grp_dim,&
                &grp_offset=grp_offset)                                         ! reuse XYZ(1)
            CHCKERR('')
            
            ! print attribute with this data item
            call print_HDF5_att(att(1),XYZ(1),'X_vec',1,.true.)
            
            ! create a grid with the topology, the geometry and the attribute
            ierr = print_HDF5_grid(grids(id),var_name,1,&
                &grid_time=0.5*pi*(id-1.0_dp)/(n_t(1)-1),grid_atts=att,&
                &reset=.true.)
            CHCKERR('')
            
            ! user output
            call writo('Finished plot for time step '//trim(i2str(id))//&
                &'/'//trim(i2str(product(n_t))))
        end do
        
        call lvl_ud(-1)
        
        ! create grid collection from individual grids and reset them
        ierr = print_HDF5_grid(time_col_grid,'time collection',2,&
            &grid_grids=grids,reset=.true.)
        CHCKERR('')
        
        ! add collection grid to HDF5 file and reset it
        call add_HDF5_item(file_info,time_col_grid,reset=.true.)
        
        ! close HDF5 file
        ierr = close_HDF5_file(file_info)
        CHCKERR('')
        
        ! deallocate
        deallocate(x_plot,y_plot,z_plot,l_plot,f_plot)
        
        call lvl_ud(-1)
    end function plot_X_vec_HDF5
end module
