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
    use grid_vars, only: grid_type
    use X_vars, only: X_2_type, modes_type, &
        &mds_X, mds_sol
    use vac_vars, only: vac_type
    use sol_vars, only: sol_type
    
    implicit none
    private
    public run_driver_sol
    
    ! global variables
#if ldebug
    logical :: debug_run_driver_sol = .false.                                   !< debug information for run_driver_sol \ldebug
    logical :: debug_interp_V = .false.                                         !< debug information for interp_v \ldebug
    integer :: n_interp_V_spline_copies                                         !< number of times a copy was done in interp_V_spline
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
            &n_procs, X_grid_style, jump_to_sol
        use grid_vars, only: n_r_sol, n_alpha
        use PB3D_ops, only: reconstruct_PB3D_grid, reconstruct_PB3D_sol
        use SLEPC_ops, only: solve_EV_system_SLEPC
        use grid_ops, only: calc_norm_range, setup_grid_sol, &
            &print_output_grid, redistribute_output_grid
        use sol_ops, only: print_output_sol
        use rich_vars, only: rich_lvl
        use rich_ops, only: calc_rich_ex
        use vac_ops, only: calc_vac_res, print_output_vac
        use grid_utilities, only: trim_grid
        use X_ops, only: setup_modes, redistribute_output_X
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
        type(grid_type) :: grid_X_rdst                                          ! redistributed perturbation grid
        type(grid_type) :: grid_X_sol                                           ! perturbation grid with solution normal part
        type(X_2_type) :: X_rdst                                                ! redistributed X
        type(X_2_type) :: X_sol                                                 ! interpolated X
        integer :: sol_limits(2)                                                ! min. and max. index of sol grid for this process
        integer :: rich_lvl_name                                                ! either the Richardson level or zero, to append to names
        real(dp), allocatable :: r_F_sol(:)                                     ! normal points in solution grid
        logical :: do_vac_ops                                                   ! whether specific calculations for vacuum are necessary
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! set up solution grid if first level
        if (rich_lvl.eq.rich_restart_lvl) then
            ! Divide solution grid under group processes, calculating the limits
            ! and the normal coordinate.
            allocate(r_F_sol(n_r_sol))
            ierr = calc_norm_range('PB3D_sol',sol_limits=sol_limits,&
                &r_F_sol=r_F_sol)
            CHCKERR('')
            
            if (rich_lvl.eq.1) then
                call writo('Set up solution grid')
                call lvl_ud(1)
                
                call writo('Calculate the grid')
                call lvl_ud(1)
                ierr = setup_grid_sol(grid_eq,grid_X,grid_sol,r_F_sol,&
                    &sol_limits)
                CHCKERR('')
                call lvl_ud(-1)
                
                call writo('Write to output file')
                call lvl_ud(1)
                ierr = print_output_grid(grid_sol,'solution','sol',&
                    &remove_previous_arrs=(jump_to_sol.and.&
                    &(X_grid_style.eq.1 .or. X_grid_style.eq.3)))
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
        
        ! initialize interpolated X
        ierr = grid_X_sol%init([1,n_alpha,grid_sol%n(3)],&
            &[grid_sol%i_min,grid_sol%i_max],grid_sol%divided)
        CHCKERR('')
        grid_X_sol%r_F = grid_sol%r_F
        grid_X_sol%r_E = grid_sol%r_E
        grid_X_sol%loc_r_F = grid_sol%loc_r_F
        grid_X_sol%loc_r_E = grid_sol%loc_r_E
        
        
        select case (X_grid_style)
            case (1,3)                                                          ! equilibrium or enriched
                ! user output
                call writo('Redistribute and interpolate the perturbation &
                    &variables to solution grid')
                call lvl_ud(1)
                
                ! redistribute grid and X variables
                ierr = redistribute_output_grid(grid_X,grid_X_rdst)
                CHCKERR('')
                ierr = redistribute_output_X(mds_X,grid_X,grid_X_rdst,X,X_rdst)
                CHCKERR('')
                
                ! initialize
                call X_sol%init(mds_sol,grid_X_sol,is_field_averaged=.true.)
                
                ! interpolate
                ierr = interp_V(mds_X,grid_X_rdst,X_rdst,mds_sol,grid_X_sol,&
                    &X_sol)
                CHCKERR('')
                
                call lvl_ud(-1)
            case (2)                                                            ! solution
                ! user output
                call writo('Copy the perturbation variables to solution grid')
                call lvl_ud(1)
                
                ! copy
                call X%copy(mds_sol,grid_X,X_sol)
                
                call lvl_ud(-1)
        end select
        
#if ldebug
        ! write integrated  field-aligned tensorial  perturbation quantities
        ! to output and plot
        if (debug_run_driver_sol) then
            ierr = print_debug_X_2(mds_sol,grid_X_sol,X_sol)
            CHCKERR('')
        end if
#endif
        
        ! clean up
        call grid_X_sol%dealloc()
        
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
        ierr = print_output_sol(grid_sol,sol,'sol',rich_lvl=rich_lvl,&
            &remove_previous_arrs=&
            &(jump_to_sol.and.(X_grid_style.eq.1 .or. X_grid_style.eq.3)))
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
        use num_vars, only: V_interp_style, rank
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
        integer :: kd                                                           ! counter
        integer :: k, m                                                         ! counters for mode numbers
        integer :: c_loc(2)                                                     ! local c for symmetric and asymmetric variables
        integer :: n_mod_tot                                                    ! total amount of modes
        integer :: kdl_i(2), kdl_o(2)                                           ! limits on normal index for a mode combination
        real(dp), pointer :: r_i_loc(:), r_o_loc(:)                             ! local r_i and r_o for a mode combination
        complex(dp), pointer :: V_i(:,:,:), V_o(:,:,:)                          ! pointers to input and output PV_i and KV_i
        logical :: calc_this(2)                                                 ! whether this combination needs to be calculated
        logical :: extrap                                                       ! whether extrapolation is used
        character(len=max_str_ln) :: err_msg                                    ! error message
#if ldebug
        integer :: km_id
        integer, allocatable :: norm_ext_i(:,:)                                 ! normal extent for input quantity mode combinations
        real(dp), allocatable :: r_loc_tot(:,:,:)                               ! r_i_loc and r_o_loc for all combinations
        complex(dp), allocatable :: V_plot(:,:)                                 ! for debug plotting of interpolated V
#endif
        
        ! initialize ierr
        ierr = 0
        
        ! test
        if (grid_i%n(2).ne.grid_o%n(2)) then
            ierr = 1
            err_msg = 'input and output grid are not compatible in the &
                &geodesic coordinate: '//trim(i2str(grid_i%n(2)))//' vs. '//&
                &trim(i2str(grid_o%n(2)))
            CHCKERR(err_msg)
        end if
        
        ! set extrapolation
        select case (V_interp_style)
            case (1)                                                            ! splines
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
        allocate(norm_ext_i(n_mod_tot,n_mod_tot))
        km_id = 0
        norm_ext_i = 0
        n_interp_V_spline_copies = 0
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
                
                ! check whether mode combination needs to be calculated
                calc_this(1) = is_necessary_X(.true.,&
                    &[mds_o%sec(k,1),mds_o%sec(m,1)])
                calc_this(2) = is_necessary_X(.false.,&
                    &[mds_o%sec(k,1),mds_o%sec(m,1)])
                
#if ldebug
                km_id = km_id + 1
                r_loc_tot(:,km_id,1) = [r_i_loc(1),r_i_loc(size(r_i_loc))]
                r_loc_tot(:,km_id,2) = [r_o_loc(1),r_o_loc(size(r_o_loc))]
                r_loc_tot(:,km_id,3) = km_id
                norm_ext_i(k,m) = size(r_i_loc)
                if (debug_interp_V) write(*,*) 'k, m', k, m, 'of', n_mod_tot, &
                        &'with calc_this = ', calc_this
#endif
                
                ! set up c_loc
                c_loc(1) = c([mds_o%sec(k,4),mds_o%sec(m,4)],.true.,n_mod_X)
                c_loc(2) = c([mds_o%sec(k,4),mds_o%sec(m,4)],.false.,n_mod_X)
                
                if (calc_this(1)) then
                    !if (k.gt.70) debug_interp_V = .true.
                    V_i => X_i%PV_0(:,:,kdl_i(1):kdl_i(2),c_loc(1))
                    V_o => X_o%PV_0(:,:,kdl_o(1):kdl_o(2),c_loc(1))
                    ierr = interp_V_spline(V_i,V_o,r_i_loc,r_o_loc,extrap)
                    CHCKERR('')
                    !debug_interp_V = .false.
                    
                    V_i => X_i%PV_2(:,:,kdl_i(1):kdl_i(2),c_loc(1))
                    V_o => X_o%PV_2(:,:,kdl_o(1):kdl_o(2),c_loc(1))
                    ierr = interp_V_spline(V_i,V_o,r_i_loc,r_o_loc,extrap)
                    CHCKERR('')
                    
                    V_i => X_i%KV_0(:,:,kdl_i(1):kdl_i(2),c_loc(1))
                    V_o => X_o%KV_0(:,:,kdl_o(1):kdl_o(2),c_loc(1))
                    ierr = interp_V_spline(V_i,V_o,r_i_loc,r_o_loc,extrap)
                    CHCKERR('')
                    
                    V_i => X_i%KV_2(:,:,kdl_i(1):kdl_i(2),c_loc(1))
                    V_o => X_o%KV_2(:,:,kdl_o(1):kdl_o(2),c_loc(1))
                    ierr = interp_V_spline(V_i,V_o,r_i_loc,r_o_loc,extrap)
                    CHCKERR('')
                end if
                
                if (calc_this(2)) then
                    V_i => X_i%PV_1(:,:,kdl_i(1):kdl_i(2),c_loc(2))
                    V_o => X_o%PV_1(:,:,kdl_o(1):kdl_o(2),c_loc(2))
                    ierr = interp_V_spline(V_i,V_o,r_i_loc,r_o_loc,extrap)
                    CHCKERR('')
                    
#if ldebug
                    if (debug_run_driver_sol) then
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
                    end if
#endif
                    
                    V_i => X_i%KV_1(:,:,kdl_i(1):kdl_i(2),c_loc(2))
                    V_o => X_o%KV_1(:,:,kdl_o(1):kdl_o(2),c_loc(2))
                    ierr = interp_V_spline(V_i,V_o,r_i_loc,r_o_loc,extrap)
                    CHCKERR('')
                end if
            end do
        end do
        
#if ldebug
        if (n_interp_V_spline_copies.gt.0) call writo('Copy was performed '//&
            &trim(i2str(n_interp_V_spline_copies))//&
            &' times in interp_V_spline',warning=.true.)
        call print_ex_2D([''],'r_i_loc_'//trim(i2str(rank)),&
            &r_loc_tot(:,1:km_id,1)/max_flux_F*2*pi,x=r_loc_tot(:,1:km_id,3),&
            &draw=.false.)
        call draw_ex([''],'r_i_loc_'//trim(i2str(rank)),km_id,1,.false.,&
            &ex_plot_style=1)
        call print_ex_2D([''],'r_o_loc_'//trim(i2str(rank)),&
            &r_loc_tot(:,1:km_id,2)/max_flux_F*2*pi,x=r_loc_tot(:,1:km_id,3),&
            &draw=.false.)
        call draw_ex([''],'r_o_loc_'//trim(i2str(rank)),km_id,1,.false.,&
            &ex_plot_style=1)
        call plot_HDF5('normal extent','norm_ext_i_'//trim(i2str(rank)),1._dp*&
            &reshape(norm_ext_i*1._dp,[n_mod_tot,n_mod_tot,1]))
#endif
        
        ! clean up
        nullify(r_i_loc,r_o_loc)
        nullify(V_i,V_o)
    contains
        !> \private Interpolation for a mode pair using splines.
        integer function interp_V_spline(V_i,V_o,r_i,r_o,extrap) result(ierr)
            
            use num_utilities, only: spline
            
            character(*), parameter :: rout_name = 'interp_V_spline'
            
            ! input / output
            complex(dp), intent(in) :: V_i(:,:,:)                               ! input tensorial metric variable for a mode combinations
            complex(dp), intent(out) :: V_o(:,:,:)                              ! output tensorial metric variable for a mode combinations
            real(dp), intent(in) :: r_i(:)                                      ! input r
            real(dp), intent(in) :: r_o(:)                                      ! output r
            logical, intent(in) :: extrap                                       ! extrapolation possible
            
            ! local variables
            integer :: id, jd, kd                                               ! counters
            character(len=max_str_ln) :: err_msg                                ! error message
            
            ! initialize ierr
            ierr = 0
            
#if ldebug
            ! test sizes
            if (size(r_i).ne.size(V_i,3)) then
                ierr = 1
                err_msg = 'r_i and V_i are not compatible'
                CHCKERR(err_msg)
            end if
            if (size(r_o).ne.size(V_o,3)) then
                ierr = 1
                err_msg = 'r_o and V_o are not compatible'
                CHCKERR(err_msg)
            end if
            if (size(V_i,1).ne.size(V_o,1) .or. size(V_i,2).ne.size(V_o,2)) then
                ierr = 1
                err_msg = 'V_i and V_o are not compatible'
                CHCKERR(err_msg)
            end if
#endif
            
            ! interpolate depending on normal size
            select case (size(r_i))
                case (1)                                                        ! copy directly
#if ldebug
                    n_interp_V_spline_copies = n_interp_V_spline_copies + 1
#endif
                    do kd = 1,size(V_o,3)
                        V_o(:,:,kd) = V_i(:,:,1)
                    end do
                case (2:3)                                                      ! linear
                    do jd = 1,size(V_i,2)
                        do id = 1,size(V_i,1)
                            ierr = spline(r_i,V_i(id,jd,:),r_o,V_o(id,jd,:),&
                                &ord=1,extrap=extrap)
                            CHCKERR('')
                        end do
                    end do
                case (4:)                                                       ! spline
                    do jd = 1,size(V_i,2)
                        do id = 1,size(V_i,1)
                            ierr = spline(r_i,V_i(id,jd,:),r_o,V_o(id,jd,:),&
                                &ord=3,extrap=extrap)
                            CHCKERR('')
                        end do
                    end do
#if ldebug
                case default
                    ierr = 1
                    CHCKERR('Need more than 0 points')
#endif
            end select
            
#if ldebug
            ! visualize
            do jd = 1,size(V_i,2)
                do id = 1,size(V_i,1)
                    if (debug_interp_V) then
                        write(*,*) '    id,jd', id, jd
                        call plot_comp('real part',r_i,rp(V_i(id,jd,:)),&
                            &r_o,rp(V_o(id,jd,:)))
                        call plot_comp('imaginary part',r_i,ip(V_i(id,jd,:)),&
                            &r_o,ip(V_o(id,jd,:)))
                    end if
                end do
            end do
#endif
        end function interp_V_spline
        
#if ldebug
        !> \private Plot to compare interpolation to original variables
        subroutine plot_comp(plot_name,x_i,y_i,x_o,y_o)
            ! input / output
            character(len=*), intent(in) :: plot_name                           ! name of plot
            real(dp), intent(in) :: x_i(:)                                      ! x input
            real(dp), intent(in) :: y_i(:)                                      ! y input
            real(dp), intent(in) :: x_o(:)                                      ! x output
            real(dp), intent(in) :: y_o(:)                                      ! y output
            
            ! local variables
            real(dp), allocatable :: x_plot(:,:)                                ! x for plot
            real(dp), allocatable :: y_plot(:,:)                                ! y for plot
            character(len=max_str_ln) :: var_names(2)                           ! variable names
            
            allocate(x_plot(max(size(x_i),size(x_o)),2))
            allocate(y_plot(max(size(x_i),size(x_o)),2))
            
            x_plot(1:size(x_i),1) = x_i
            x_plot(size(x_i):size(x_plot,1),1) = x_i(size(x_i))
            x_plot(1:size(x_o),2) = x_o
            x_plot(size(x_o):size(x_plot,1),2) = x_o(size(x_o))
            
            y_plot(1:size(x_i),1) = y_i(:)
            y_plot(size(x_i):size(y_plot,1),1) = y_i(size(x_i))
            y_plot(1:size(x_o),2) = y_o(:)
            y_plot(size(x_o):size(y_plot,1),2) = y_o(size(x_o))
            
            var_names(1) = 'input '//plot_name
            var_names(2) = 'output '//plot_name
            call print_ex_2D(var_names,'',y_plot,x=x_plot)
        end subroutine plot_comp
#endif
    end function interp_V
end module driver_sol
