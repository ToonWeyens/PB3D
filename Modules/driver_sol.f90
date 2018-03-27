!------------------------------------------------------------------------------!
!> Driver of the solution part of PB3D.
!------------------------------------------------------------------------------!
module driver_sol
#include <PB3D_macros.h>
#include <wrappers.h>
    use str_utilities
    use output_ops
    use messages
    use num_vars, only: dp, pi, max_str_ln, iu
    use grid_vars, only: grid_type, disc_type
    use X_vars, only: X_2_type, modes_type, &
        &mds_X, mds_sol
    use vac_vars, only: vac_type
    use sol_vars, only: sol_type
    
    implicit none
    private
    public run_driver_sol
    
    ! global variables
#if ldebug
    logical :: debug_run_driver_sol = .true.                                   !< debug information for run_driver_sol \ldebug
#endif
    
contains
    !> Main driver of PB3D solution part.
    !!
    !!  - sets up:
    !!      * \c grid_sol (only first Richardson level)
    !!      * \c sol
    !!  - writes to HDF5:
    !!      * \c grid_sol (only first Richardson level)
    !!      * \c sol
    !!  - deallocates:
    !!      * sol before setting up (but after guess)
    !!
    !! \return ierr
    integer function run_driver_sol(grid_eq,grid_X,grid_sol,X,vac,sol) &
        &result(ierr)
        
        use num_vars, only: EV_style, eq_style, rich_restart_lvl, rank, &
            &n_procs, X_grid_style, norm_disc_prec_X
        use grid_vars, only: n_r_sol, n_alpha
        use PB3D_ops, only: reconstruct_PB3D_grid, reconstruct_PB3D_sol
        use SLEPC_ops, only: solve_EV_system_SLEPC
        use grid_ops, only: calc_norm_range, setup_grid_sol, &
            &print_output_grid
        use sol_ops, only: print_output_sol
        use rich_vars, only: rich_lvl
        use rich_ops, only: calc_rich_ex
        use vac_ops, only: calc_vac_res, print_output_vac
        use MPI_utilities, only: get_ser_var
        use grid_utilities, only: trim_grid, setup_interp_data, apply_disc
        use X_ops, only: setup_modes
        use X_vars, only: min_m_X, min_n_X
#if ldebug
        use X_ops, only: print_debug_X_2
#endif
        
        character(*), parameter :: rout_name = 'run_driver_sol'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  !< equilibrium grid
        type(grid_type), intent(in), target :: grid_X                           !< perturbation grid
        type(grid_type), intent(inout) :: grid_sol                              !< solution grid
        type(X_2_type), intent(in) :: X                                         !< integrated tensorial perturbation variables
        type(vac_type), intent(inout) :: vac                                    !< vacuum variables
        type(sol_type), intent(inout) :: sol                                    !< solution variables
        
        ! local variables
        type(grid_type) :: grid_sol_trim                                        ! trimmed solution grid
        type(grid_type) :: grid_X_trim                                          ! trimmed perturbation grid
        type(grid_type) :: grid_X_ser                                           ! serial perturbation grid
        type(grid_type) :: grid_X_sol                                           ! perturbation grid with solution normal part
        type(X_2_type) :: X_ser                                                 ! serial X
        type(X_2_type) :: X_sol                                                 ! interpolated X
        integer :: pmone                                                        ! plus (increasing r_F) or minus (decreasing r_F) one
        integer :: min_nm_X                                                     ! minimal n (tor. flux) or m (pol. flux)
        integer :: kd, ld                                                       ! counters
        integer :: n_mod_tot                                                    ! total number of modes
        integer :: norm_range_X(2)                                              ! normal range in X grid of total mode
        integer :: norm_range_sol(2)                                            ! normal range in sol grid of total mode
        integer :: n_X(3)                                                       ! n of grid_X_trim
        integer :: n_X_loc                                                      ! local size in grid_X_trim
        integer :: sol_limits(2)                                                ! min. and max. index of sol grid for this process
        integer :: norm_id(2)                                                   ! untrimmed normal indices for trimmed solution grid
        integer :: rich_lvl_name                                                ! either the Richardson level or zero, to append to names
        real(dp), allocatable :: r_F_sol(:)                                     ! normal points in solution grid
        complex(dp), allocatable :: V_X(:,:,:)                                  ! PV or KV for total mode
        complex(dp), allocatable :: ser_var_loc(:)                              ! local serial variable
        logical :: do_vac_ops                                                   ! whether specific calculations for vacuum are necessary
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! set up solution grid if first level
        if (rich_lvl.eq.rich_restart_lvl) then
            ! Divide solution grid under group processes, calculating the limits
            ! and the normal coordinate.
            allocate(r_F_sol(n_r_sol))
            ierr = calc_norm_range(sol_limits=sol_limits,r_F_sol=r_F_sol)
            CHCKERR('')
            
            if (rich_lvl.eq.1) then
                call writo('Set up solution grid')
                call lvl_ud(1)
                
                call writo('Calculate the grid')
                call lvl_ud(1)
                ierr = setup_grid_sol(grid_X,grid_sol,r_F_sol,sol_limits)
                CHCKERR('')
                call lvl_ud(-1)
                
                call writo('Write to output file')
                call lvl_ud(1)
                ierr = print_output_grid(grid_sol,'solution','sol')
                CHCKERR('')
                call lvl_ud(-1)
                
                ! set up modes
                ierr = setup_modes(mds_sol,grid_eq,grid_sol,plot_nm=.false.)
                CHCKERR('')
                
                call lvl_ud(-1)
            else
                ! restore solution grid and trim it
                ierr = reconstruct_PB3D_grid(grid_sol,'sol',&
                    &grid_limits=sol_limits)
                CHCKERR('')
                ierr = trim_grid(grid_sol,grid_sol_trim)
                CHCKERR('')
                
                ! reconstruct solution on trimmed grid
                ierr = reconstruct_PB3D_sol(mds_sol,grid_sol_trim,sol,'sol',&
                    &rich_lvl=rich_lvl-1)
                CHCKERR('')
                
                ! clean up
                call grid_sol_trim%dealloc()
            end if
            
            deallocate(r_F_sol)
        end if
        
        ! set up  whether Richardson level  has to be  appended to the  name and
        ! whether to do vacuum operations
        select case (eq_style) 
            case (1)                                                            ! VMEC
                rich_lvl_name = rich_lvl                                        ! append richardson level
                do_vac_ops = .true.
            case (2)                                                            ! HELENA
                rich_lvl_name = 0                                               ! do not append name
                if (rich_lvl.eq.rich_restart_lvl) then
                    do_vac_ops = .true.
                else
                    do_vac_ops = .false.
                end if
        end select
        write(*,*) '¡¡¡¡¡ NO VACUUM !!!!!'
        do_vac_ops = .false.
        
        if (do_vac_ops) then
            ! calculate vacuum
            ierr = calc_vac_res(mds_sol,vac)
            CHCKERR('')
            
            call writo('Write to output file')
            call lvl_ud(1)
            CHCKERR('')
            if (rank.eq.n_procs-1) then
                ierr = print_output_vac(vac,'vac',rich_lvl=rich_lvl_name)
                CHCKERR('')
            end if
            call lvl_ud(-1)
        end if
        
        select case (X_grid_style)
            case (1)                                                            ! equilibrium
                ! user output
                call writo('Interpolate the perturbation variables to &
                    &solution grid')
                call lvl_ud(1)
                
                ! initialize interpolated X
                ierr = grid_X_sol%init([1,n_alpha,grid_sol%n(3)],&
                    &[grid_sol%i_min,grid_sol%i_max],grid_sol%divided)
                CHCKERR('')
                grid_X_sol%r_F = grid_sol%r_F
                grid_X_sol%r_E = grid_sol%r_E
                grid_X_sol%loc_r_F = grid_sol%loc_r_F
                grid_X_sol%loc_r_E = grid_sol%loc_r_E
                call X_sol%init(mds_sol,grid_X_sol,is_field_averaged=.true.)
                
                ! interpolate
                ierr = interp_V(mds_X,grid_X,X,mds_sol,grid_X_sol,X_sol)
                CHCKERR('')
                
#if ldebug
                ! write  integrated  and  interpolated  field-aligned  tensorial
                ! perturbation quantities to output and plot
                if (debug_run_driver_sol) then
                    ierr = print_debug_X_2(mds_sol,grid_X_sol,X_sol)
                    CHCKERR('')
                end if
#endif
                
                ! clean up
                call grid_X_sol%dealloc()
                
                call lvl_ud(-1)
            case (2)                                                            ! solution
                ! user output
                call writo('Copy the perturbation variables to solution grid')
                call lvl_ud(1)
                
                ! copy
                call X%copy(mds_sol,grid_X,X_sol)
                
                call lvl_ud(-1)
        end select
        
        ! solve the system
        call writo('Solving the system')
        call lvl_ud(1)
        select case (EV_style)
            case(1)                                                             ! SLEPC solver for EV problem
                ! solve the system
                ierr = solve_EV_system_SLEPC(mds_sol,grid_sol,X_sol,vac,sol)
                CHCKERR('')
            case default
                err_msg = 'No EV solver style associated with '//&
                    &trim(i2str(EV_style))
                ierr = 1
                CHCKERR(err_msg)
        end select
        call lvl_ud(-1)
        call writo('System solved')
        
        ! write solution variables to output
        ierr = print_output_sol(grid_sol,sol,'sol',rich_lvl=rich_lvl)
        CHCKERR('')
        
        ! calculate Richardson extrapolation factors if necessary
        call calc_rich_ex(sol%val)
        
        ! clean up
        call X_sol%dealloc()
    end function run_driver_sol
    
    !> \public  Interpolate  tensorial  perturbation  quantities  in  the  third
    !! dimension.
    !!
    !! The input grid should not be divided, whereas the output grid can be.
    !!
    !! The procedure  considers all  possible mode  number combinations.  For \c
    !! X_style 2  (fast), each  secondary mode  only lives  in a  certain normal
    !! range  of the  plasma. Therefore,  each secondary  mode pair  also has  a
    !! limited normal range, given by the overlap of the ranges of the members.
    !!
    !! The interpolated mode number combinations  have a normal range that might
    !! slightly differ from the input ranges, in which case extrapolation can be
    !! done if the method allows for it.
    !!
    !! \note If the input grid is too  coarse, but the interpolated grid is not,
    !! it is in  theory possible that there are mode  numbers and therefore mode
    !! number pairs  that exist  in the  input grid,  though they  do so  in the
    !! interpolated  grid.  In  this  case,  they  can  not  be  calculated  and
    !! are  currently  set  to  zero. However,  the  current  implementation  of
    !! setup_modes()  produces an  error if  less than  the full  range of  mode
    !! numbers is accurately present.
    integer function interp_V(mds_i,grid_i,X_i,mds_o,grid_o,X_o) result(ierr)
        use X_utilities, only: is_necessary_X
        use X_vars, only: n_mod_X
        use num_vars, only: V_interp_style, norm_disc_prec_sol
        use grid_utilities, only: setup_interp_data
        use eq_vars, only: max_flux_F
        use num_utilities, only: c
        
        character(*), parameter :: rout_name = 'interp_V'
        
        ! input / output
        type(modes_type), intent(in) :: mds_i                                   !< general modes variables for input
        type(grid_type), intent(in) :: grid_i                                   !< grid at which \c X_i is tabulated
        type(X_2_type), intent(in), target :: X_i                               !< tensorial perturbation variable on input grid
        type(modes_type), intent(in) :: mds_o                                   !< general modes variables for output
        type(grid_type), intent(in) :: grid_o                                   !< grid at which \c X_o is interpolated
        type(X_2_type), intent(inout), target :: X_o                            !< interpolated tensorial perturbation variable
        
        ! local variables
        type(disc_type) :: norm_interp_data                                     ! data for normal interpolation X->sol
        integer :: id, jd, kd, md                                               ! counter
        integer :: k, m                                                         ! counters for mode numbers
        integer :: c_loc(2)                                                     ! local c for symmetric and asymmetric variables
        integer :: n_mod_tot                                                    ! total amount of modes
        integer :: norm_disc_prec_loc                                           ! local precision
        integer :: kdl_i(2), kdl_o(2)                                           ! limits on normal index for a mode combination
        real(dp), pointer :: r_i_loc(:), r_o_loc(:)                             ! local r_i and r_o for a mode combination
        real(dp), allocatable :: spline_knots(:,:)                              ! knots of spline
        real(dp), allocatable :: spline_coeff(:,:)                              ! coefficients of spline
        complex(dp), pointer :: V_i(:,:,:), V_o(:,:,:)                          ! pointers to input and output PV_i and KV_i
        logical :: calc_this(2)                                                 ! whether this combination needs to be calculated
        logical :: extrap                                                       ! whether extrapolation is used
        character(len=max_str_ln) :: err_msg                                    ! error message
#if ldebug
        integer :: km_id
        real(dp), allocatable :: r_loc_tot(:,:,:)                               ! r_i_loc and r_o_loc for all combinations
        complex(dp), allocatable :: V_plot(:,:)                                 ! for debug plotting of interpolated V
#endif
        
        ! initialize ierr
        ierr = 0
        
        ! test
        if (grid_i%n(2).ne.grid_o%n(2)) then
            ierr = 1
            err_msg = 'input and output grid are not compatible in the geodesic &
                &coordinate: '//trim(i2str(grid_i%n(2)))//' vs. '//&
                &trim(i2str(grid_o%n(2)))
            CHCKERR(err_msg)
        end if
        
        ! set extrapolation
        select case (V_interp_style)
            case (1)                                                            ! finite differences
                !extrap = .false.
                extrap = .true.
            case (2)                                                            ! splines
                extrap = .true.
        end select
        
        ! initialize to zero
        X_o%PV_0 = 0._dp
        X_o%PV_1 = 0._dp
        X_o%PV_2 = 0._dp
        X_o%KV_0 = 0._dp
        X_o%KV_1 = 0._dp
        X_o%KV_2 = 0._dp
        
        ! select all mode number combinations of interpolated grid
        n_mod_tot = size(mds_o%sec,1)
#if ldebug
        allocate(r_loc_tot(2,n_mod_tot**2,3))
        km_id = 0
#endif
        do m = 1,n_mod_tot
            do k = 1,n_mod_tot
                ! test whether input and output grid are consistent
                if (mds_i%sec(k,4).ne.mds_o%sec(k,4) .or. &
                    & mds_i%sec(m,4).ne.mds_o%sec(m,4)) then
                    ierr = 1
                    err_msg = 'For ('//trim(i2str(k))//','//trim(i2str(m))//&
                        &'), no consistency in grid_i and grid_o'
                    CHCKERR(err_msg)
                end if
                
                ! set input and output grid normal limits for mode pair (k,m)
                kdl_i(1) = max(mds_i%sec(k,2),mds_i%sec(m,2))
                kdl_i(2) = min(mds_i%sec(k,3),mds_i%sec(m,3))
                kdl_o(1) = max(mds_o%sec(k,2),mds_o%sec(m,2))
                kdl_o(2) = min(mds_o%sec(k,3),mds_o%sec(m,3))
                
                ! limit to grid ranges
                kdl_i(1) = max(kdl_i(1),grid_i%i_min)
                kdl_i(2) = min(kdl_i(2),grid_i%i_max)
                kdl_o(1) = max(kdl_o(1),grid_o%i_min)
                kdl_o(2) = min(kdl_o(2),grid_o%i_max)
                
                ! limit output range to input range if not extrapolating
                if (.not.extrap) then
                    do kd = kdl_o(1),kdl_o(2)
                        if (grid_o%r_F(kd).ge.grid_i%r_F(kdl_i(1))) exit
                    end do
                    kdl_o(1) = kd
                    
                    do kd = kdl_o(2),kdl_o(1),-1
                        if (grid_o%r_F(kd).le.grid_i%r_F(kdl_i(2))) exit
                    end do
                    kdl_o(2) = kd
                end if
                
                ! cycle if mode combination does not exist anywhere
                if (kdl_i(1).gt.kdl_i(2) .or. kdl_o(1).gt.kdl_o(2)) cycle
                
                ! convert limits to local
                kdl_i = kdl_i - grid_i%i_min + 1
                kdl_o = kdl_o - grid_o%i_min + 1
                
                ! get input ranges for this mode number
                r_i_loc => grid_i%loc_r_F(kdl_i(1):kdl_i(2))
                r_o_loc => grid_o%loc_r_F(kdl_o(1):kdl_o(2))
                
#if ldebug
                km_id = km_id + 1
                r_loc_tot(:,km_id,1) = [r_i_loc(1),r_i_loc(size(r_i_loc))]
                r_loc_tot(:,km_id,2) = [r_o_loc(1),r_o_loc(size(r_o_loc))]
                r_loc_tot(:,km_id,3) = km_id
#endif
                
                ! check whether mode combination needs to be calculated
                calc_this(1) = is_necessary_X(.true.,&
                    &[mds_o%sec(k,1),mds_o%sec(m,1)])
                calc_this(2) = is_necessary_X(.false.,&
                    &[mds_o%sec(k,1),mds_o%sec(m,1)])
                
                ! set up c_loc
                c_loc(1) = c([mds_o%sec(k,4),mds_o%sec(m,4)],.true.,n_mod_X)
                c_loc(2) = c([mds_o%sec(k,4),mds_o%sec(m,4)],.false.,n_mod_X)
                
                ! set precision
                norm_disc_prec_loc = min(kdl_i(2)-kdl_i(1),&
                    &norm_disc_prec_sol)
                
                ! prepare interpolation
                select case (V_interp_style)
                    case (1)                                                    ! finite differences
                        ! can only interpolate with at least precision one
                        !if (norm_disc_prec_loc.ge.1) then
                            ! set up normal interpolation factors
                            ierr = setup_interp_data(r_i_loc,r_o_loc,&
                                &norm_interp_data,norm_disc_prec_loc,&
                                &extrap=.true.)
                            CHCKERR('')
                        !else
                            !write(*,*) 'ONLY ONE POINT'
                            !read(*,*)
                        !end if
                    case (2)                                                    ! splines
                        allocate(spline_coeff(size(r_i_loc),2))
                        allocate(spline_knots(size(r_i_loc)+norm_disc_prec_loc,&
                            &2))
                end select
                
                if (calc_this(1)) then
                    V_i => X_i%PV_0(:,:,kdl_i(1):kdl_i(2),c_loc(1))
                    V_o => X_o%PV_0(:,:,kdl_o(1):kdl_o(2),c_loc(1))
                    ierr = interp_V_loc(V_i,V_o,norm_interp_data,&
                        &spline_knots,spline_coeff,norm_disc_prec_loc,extrap,&
                        &V_interp_style)
                    CHCKERR('')
                    
                    V_i => X_i%PV_2(:,:,kdl_i(1):kdl_i(2),c_loc(1))
                    V_o => X_o%PV_2(:,:,kdl_o(1):kdl_o(2),c_loc(1))
                    ierr = interp_V_loc(V_i,V_o,norm_interp_data,&
                        &spline_knots,spline_coeff,norm_disc_prec_loc,extrap,&
                        &V_interp_style)
                    CHCKERR('')
                    
                    V_i => X_i%KV_0(:,:,kdl_i(1):kdl_i(2),c_loc(1))
                    V_o => X_o%KV_0(:,:,kdl_o(1):kdl_o(2),c_loc(1))
                    ierr = interp_V_loc(V_i,V_o,norm_interp_data,&
                        &spline_knots,spline_coeff,norm_disc_prec_loc,extrap,&
                        &V_interp_style)
                    CHCKERR('')
                    
                    V_i => X_i%KV_2(:,:,kdl_i(1):kdl_i(2),c_loc(1))
                    V_o => X_o%KV_2(:,:,kdl_o(1):kdl_o(2),c_loc(1))
                    ierr = interp_V_loc(V_i,V_o,norm_interp_data,&
                        &spline_knots,spline_coeff,norm_disc_prec_loc,extrap,&
                        &V_interp_style)
                    CHCKERR('')
                end if
                
                if (calc_this(2)) then
                    V_i => X_i%PV_1(:,:,kdl_i(1):kdl_i(2),c_loc(2))
                    V_o => X_o%PV_1(:,:,kdl_o(1):kdl_o(2),c_loc(2))
                    ierr = interp_V_loc(V_i,V_o,norm_interp_data,&
                        &spline_knots,spline_coeff,norm_disc_prec_loc,extrap,&
                        &V_interp_style)
                    CHCKERR('')
                    
#if ldebug
                    !if (debug_run_driver_sol) then
                        !call writo('For [k,m] = ['//&
                            !&trim(i2str(mds_o%sec(k,1)))//','//&
                            !&trim(i2str(mds_o%sec(m,1)))//']:')
                        !allocate(V_plot(max(size(V_i,3),size(V_o,3)),4))
                        !V_plot(1:size(V_i,3),1) = V_i(1,1,:)
                        !V_plot(1:size(V_o,3),2) = V_o(1,1,:)
                        !V_plot(1:size(V_i,3),3) = r_i_loc/max_flux_F*2*pi
                        !V_plot(1:size(V_o,3),4) = r_o_loc/max_flux_F*2*pi
                        !V_plot(size(V_i,3)+1:size(V_plot,1),1) = &
                            !&V_plot(size(V_i,3),1)
                        !V_plot(size(V_o,3)+1:size(V_plot,1),2) = &
                            !&V_plot(size(V_o,3),2)
                        !V_plot(size(V_i,3)+1:size(V_plot,1),3) = &
                            !&V_plot(size(V_i,3),3)
                        !V_plot(size(V_o,3)+1:size(V_plot,1),4) = &
                            !&V_plot(size(V_o,3),4)
                        !call print_ex_2D(['V_i','V_o'],'RE_V',&
                            !&rp(V_plot(:,1:2)),x=rp(V_plot(:,3:4)))
                        !!call print_ex_2D(['V_i','V_o'],'IM_V',&
                            !!&ip(V_plot(:,1:2)),x=rp(V_plot(:,3:4)))
                        !deallocate(V_plot)
                    !end if
#endif
                    
                    V_i => X_i%KV_1(:,:,kdl_i(1):kdl_i(2),c_loc(2))
                    V_o => X_o%KV_1(:,:,kdl_o(1):kdl_o(2),c_loc(2))
                    ierr = interp_V_loc(V_i,V_o,norm_interp_data,&
                        &spline_knots,spline_coeff,norm_disc_prec_loc,extrap,&
                        &V_interp_style)
                    CHCKERR('')
                end if
                
                ! clean up
                select case (V_interp_style)
                    case (1)                                                    ! finite differences
                        if (norm_disc_prec_loc.ge.1) then
                            call norm_interp_data%dealloc()
                        end if
                    case (2)                                                    ! splines
                        deallocate(spline_knots)
                        deallocate(spline_coeff)
                end select
            end do
        end do
        
#if ldebug
        call print_ex_2D([''],'r_i_loc',&
            &r_loc_tot(:,1:km_id,1)/max_flux_F*2*pi,x=r_loc_tot(:,1:km_id,3),&
            &draw=.false.)
        call draw_ex([''],'r_i_loc',km_id,1,.false.)
        call print_ex_2D([''],'r_o_loc',&
            &r_loc_tot(:,1:km_id,2)/max_flux_F*2*pi,x=r_loc_tot(:,1:km_id,3),&
            &draw=.false.)
        call draw_ex([''],'r_o_loc',km_id,1,.false.)
#endif
        
        ! clean up
        nullify(r_i_loc,r_o_loc)
        nullify(V_i,V_o)
    contains
        !> \private Interpolation for a mode pair
        !!
        !! The  specific variables  for  a style,  such  as \c  norm_interp_data
        !! always have to be passed, also when using another style, but they can
        !! be empty.
        integer function interp_V_loc(V_i,V_o,norm_interp_data,spline_knots,&
            &spline_coeff,norm_disc_prec,extrap,style) result(ierr)
            
            use grid_utilities, only: apply_disc
            use bspline_sub_module, only: db1ink, db1val, get_status_message
            
            character(*), parameter :: rout_name = 'interp_V_loc'
            
            ! input / output
            complex(dp), intent(in) :: V_i(:,:,:)                               ! input tensorial metric variable for a mode combinations
            complex(dp), intent(inout) :: V_o(:,:,:)                            ! output tensorial metric variable for a mode combinations
            type(disc_type) :: norm_interp_data                                 ! normal interpolation data
            real(dp), intent(inout) :: spline_knots(:,:)                        ! knots of spline
            real(dp), intent(inout) :: spline_coeff(:,:)                        ! coefficients of spline
            integer, intent(in) :: norm_disc_prec                               ! normal discretization
            logical, intent(in) :: extrap                                       ! extrapolation possible
            integer, intent(in) :: style                                        ! style for V_interp
            
            ! local variables
            integer :: id, jd, kd, md                                           ! counters
            integer :: spline_init(2)                                           ! spline initialization parameter
            real(dp) :: V_loc(2)                                                ! local PV or KV value
            
            ! initialize ierr
            ierr = 0
            
            if (norm_disc_prec.ge.1) then                                       ! interpolate
                select case (style)
                    case (1)                                                    ! finite differences
                        ierr = apply_disc(V_i,norm_interp_data,V_o,3)
                        CHCKERR('')
                    case (2)                                                    ! splines
                    do jd = 1,size(V_o,2)
                        do id = 1,size(V_o,1)
                            call db1ink(r_i_loc,size(r_i_loc),rp(V_i(id,jd,:)),&
                                &norm_disc_prec,0,spline_knots(:,1),&
                                &spline_coeff(:,1),ierr)
                            err_msg = get_status_message(ierr)
                            CHCKERR(err_msg)
                            call db1ink(r_i_loc,size(r_i_loc),ip(V_i(id,jd,:)),&
                                &norm_disc_prec,0,spline_knots(:,2),&
                                &spline_coeff(:,2),ierr)
                            err_msg = get_status_message(ierr)
                            CHCKERR(err_msg)
                            spline_init = 1
                            do kd = 1,size(r_o_loc)
                                do md = 1,2
                                    call db1val(r_o_loc(kd),0,&
                                        &spline_knots(:,md),size(r_i_loc),&
                                        &norm_disc_prec,spline_coeff(:,md),&
                                        &V_loc(md),ierr,spline_init(md),&
                                        &extrap=extrap)
                                    err_msg = get_status_message(ierr)
                                    CHCKERR(err_msg)
                                end do
                                V_o(id,jd,kd) = V_loc(1) + iu*V_loc(2)
                            end do
                        end do
                    end do
                end select
            else                                                                ! just copy
                do kd = 1,kdl_o(2)-kdl_o(1)+1
                    V_o(:,:,kd) = V_i(:,:,1)
                end do
            end if
        end function interp_V_loc
    end function interp_V
end module driver_sol
