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
    
    ! plots an eigenfunction in either 1D (normal coord.) or 2D (normal coord. &
    ! par coord.), depending on whether the variable par_range_F is provided
    ! The  individual modes are  plotted at t  = 0 and  t = (pi/2)/omega  and an
    ! animated plot is produced as well in a .gif output file
    ! [MPI] Collective call
    integer function plot_X_vec(X_vec,X_val,X_id,par_range_F) result(ierr)
        use num_vars, only: grp_rank, min_r_X, max_r_X, max_n_plots, &
            &grp_n_procs, use_pol_flux
        use output_ops, only: draw_GP, draw_GP_animated
        use X_vars, only: n_r_X, size_X, n_X, m_X, grp_r_X, grp_min_r_X
        use VMEC_vars, only: n_r_eq, VMEC_use_pol_flux
        use MPI_ops, only: get_ser_X_vec, get_ghost_X_vec
        use eq_vars, only: q_saf_FD, rot_t_FD, A_0, grp_n_r_eq, grp_min_r_eq, &
            &max_flux, max_flux_VMEC, flux_p_FD, flux_t_FD
        use utilities, only: dis2con, con2dis, interp_fun_1D
        
        character(*), parameter :: rout_name = 'plot_X_vec'
        
        ! input / output
        complex(dp), intent(in) :: X_vec(:,:)                                   ! MPI Eigenvector
        complex(dp), intent(in) :: X_val                                        ! Eigenvalue
        integer, intent(in) :: X_id                                             ! nr. of Eigenvalue (for output name)
        real(dp), intent(in), optional :: par_range_F(2)                        ! the parallel range [2pi]
        
        ! local variables
        real(dp), allocatable :: x_plot(:,:,:,:)                                ! x values of plot, parallel version
        real(dp), allocatable :: y_plot(:,:,:,:)                                ! y values of plot, parallel version
        real(dp), allocatable :: z_plot(:,:,:,:)                                ! z values of plot, parallel version
        real(dp), allocatable :: z_magn_plot(:,:,:)                             ! z values of plot of magnitude, parallel version
        integer :: id, jd, kd, ld                                               ! counter
        complex(dp) :: omega                                                    ! eigenvalue
        integer :: n_t                                                          ! nr. of time steps in quarter period
        integer :: n_n_t                                                        ! how many quarter periods
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
        integer :: grp_n_norm                                                   ! nr. of normal points in plot for this rank
        complex(dp) :: fac_i_interp(1)                                          ! complex copy of fac_n_interp or fac_m_interp
        real(dp), allocatable :: ranges(:,:)                                    ! ranges for animated plot
        
        ! initialize ierr
        ierr = 0
        
        call lvl_ud(1)
        
        ! 1. initialize some quantities
        
        ! set up fac_n and fac_m
        allocate(fac_n(grp_n_r_eq),fac_m(grp_n_r_eq))
        if (use_pol_flux) then
            fac_n = q_saf_FD(:,0) * A_0                                         ! n_0 = 1/A_0
            fac_m = 1.0_dp
        else
            fac_n = 1.0_dp
            fac_m = rot_t_FD(:,0)
        end if
        
        ! allocate ranges if group master
        if (grp_rank.eq.0) then
            if (present(par_range_F)) then
                allocate(ranges(3,2))
            else
                allocate(ranges(2,2))
            end if
        end if
        
        ! if  omega  is  predominantly  real,  the  Eigenfunction  explodes,
        ! so  limit  n_n_t to  1.  If  it  is predominantly  imaginary,  the
        ! Eigenfunction oscillates, so choose n_n_t = 8 for 2 whole periods
        omega = sqrt(X_val)
        if (abs(realpart(omega)/imagpart(omega)).gt.1) then
            n_n_t = 1                                                           ! 1 quarter period
            n_t = 20                                                            ! 20 points per quarter period
        else
            n_n_t = 8                                                           ! 8 quarter periods
            n_t = 5                                                             ! 5 points per quarter period
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
        
        ! 2. calculate the x and z parts of the plot in parallel
        
        ! calculate starting and ending normal index of array
        call calc_r_plot(r_plot,grp_n_norm)
        
        ! allocate variables
        allocate(x_plot(grp_n_norm,n_par_F,size_X,n_n_t*n_t))
        allocate(y_plot(grp_n_norm,n_par_F,size_X,n_n_t*n_t))
        allocate(z_plot(grp_n_norm,n_par_F,size_X,n_n_t*n_t))
        allocate(z_magn_plot(grp_n_norm,n_par_F,n_n_t*n_t))
        allocate(X_vec_interp(size_X))
        allocate(first_X_vec_interp(size_X))
        
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
        do id = 1,n_n_t*n_t
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
                    X_vec_interp = X_vec(:,floor(kd_loc_i)) + &
                        &(kd_loc_i-floor(kd_loc_i)) * &
                        &(X_vec(:,ceiling(kd_loc_i))-X_vec(:,floor(kd_loc_i)))
                    
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
                            &exp(omega/abs(omega)*(id-1.0_dp)/(n_t*n_n_t-1)&    ! time dependency
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
                plot_title = 'EV '//trim(i2str(X_id))//&
                    &' - mode '//trim(i2str(ld))//' of '//trim(i2str(size_X))
                file_name = 'mode_'//trim(i2str(X_id))//'_'//&
                    &trim(i2str(ld))//'_'//trim(i2str(grp_rank))
                if (present(par_range_F)) then
                    call print_GP_3D(trim(plot_title),trim(file_name),&
                        &z_plot(:,:,ld,:),y=y_plot(:,:,ld,:),&
                        &x=x_plot(:,:,ld,:),draw=.false.)
                else
                    call print_GP_2D(trim(plot_title),trim(file_name),&
                        &z_plot(:,1,ld,:),x=x_plot(:,1,ld,:),&
                        &draw=.false.)
                end if
                    
                ! plot by group master
                if (grp_rank.eq.0) then
                    ! set up file names in array
                    allocate(file_names(grp_n_procs))
                    do id = 1,grp_n_procs
                        file_names(id) = 'mode_'//trim(i2str(X_id))//&
                            &'_'//trim(i2str(ld))//'_'//trim(i2str(id-1))
                    end do
                    
                    ! draw animation
                    call draw_GP_animated(trim(plot_title),file_names,n_n_t*n_t,&
                        &.false.)
                    
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
        plot_title = 'EV '//trim(i2str(X_id))//' - global mode'
        file_name = 'global_mode_'//trim(i2str(X_id))//'_'//&
            &trim(i2str(grp_rank))
        if (present(par_range_F)) then
            call print_GP_3D(trim(plot_title),trim(file_name),&
                &z_magn_plot(:,:,:),y=y_plot(:,:,1,:),&
                &x=x_plot(:,:,1,:),draw=.false.)
        else
            call print_GP_2D(trim(plot_title),trim(file_name),&
                &z_magn_plot(:,1,:),x=x_plot(:,1,1,:),&
                &draw=.false.)
        end if
        
        ! plot by group master
        if (grp_rank.eq.0) then
            ! set up file names in array
            allocate(file_names(grp_n_procs))
            do id = 1,grp_n_procs
                file_names(id) = 'global_mode_'//trim(i2str(X_id))&
                    &//'_'//trim(i2str(id-1))
            end do
            
            ! draw animation
            call draw_GP_animated(trim(plot_title),file_names,n_n_t*n_t,&
                &.false.)
            
            ! deallocate
            deallocate(file_names)
        end if
        
        call lvl_ud(-1)
        
    contains
        ! calculates the array of normal  points where to interpolate to produce
        ! the plot
        subroutine calc_r_plot(r_plot,grp_n_norm)
            ! input / output
            real(dp), intent(inout), allocatable :: r_plot(:)                   ! local normal values at which to interpolate for the plot
            integer, intent(inout) :: grp_n_norm                                ! nr. of normal points for this rank
            
            ! local variables
            integer :: start_i_norm, end_i_norm                                 ! starting and ending index of normal points for this rank
            integer :: max_n_norm = 50                                          ! max nr of normal points to plot
            integer :: n_norm                                                   ! nr. of normal points
            real(dp) :: inc_norm                                                ! normal increment for plots
            integer :: kd                                                       ! counter
            
            ! set  up how  many  normal  points for  this  rank  n_norm and  the
            ! increment between normal points inc_norm
            ! (taking all the normal points could be way too many)
            n_norm = min(n_r_X,max_n_norm)
            inc_norm = (max_r_X-min_r_X)/(n_norm-1)
            
            ! calculate the  starting and ending index,  first determining which
            ! is the range where normal interpolation can take place
            start_i_norm = ceiling((grp_r_X(1)-min_r_X)/inc_norm + 1)
            end_i_norm = floor((grp_r_X(size(grp_r_X))-min_r_X)/inc_norm + 1)   ! size(grp_r_X) can be grp_n_r_X or grp_n_r_X + 1 (ghost region in X grid)
            if (grp_rank+1.lt.grp_n_procs) end_i_norm = end_i_norm + 1
            
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
    end function plot_X_vec
end module
