!------------------------------------------------------------------------------!
!> Operations that have to do with the grids and different coordinate systems.
!------------------------------------------------------------------------------!
module grid_ops
#include <PB3D_macros.h>
    use str_utilities
    use output_ops
    use messages
    use num_vars, only: dp, pi, max_str_ln
    use grid_vars, only: grid_type, disc_type
    use eq_vars, only: eq_1_type

    implicit none
    private
    public calc_norm_range, calc_ang_grid_eq_B, magn_grid_plot, setup_grid_eq, &
        &setup_grid_eq_B, setup_grid_X, setup_grid_sol, print_output_grid, &
        &redistribute_output_grid
#if ldebug
    public debug_calc_ang_grid_eq_B
#endif
    
    ! global variables
#if ldebug
    !> \ldebug
    logical :: debug_calc_ang_grid_eq_B = .false.                               !< plot debug information for calc_ang_grid_eq_b()
#endif

contains
    !> Calculates normal range  for the input grid, the  equilibrium grid and/or
    !! the solution grid.
    !!
    !! General workings, depending on \c X_grid_style:
    !!
    !! | 1 (equilibrium)     | 2 (solution )       | 3 (enriched)           |
    !! | ------------------- | ------------------- | ---------------------- |
    !! | calc_eq()           | calc_eq()           | calc_eq()              |
    !! |                     | redistribute to sol | add points to optimize |
    !! |                     | interpolate to sol  | interpolate to X       |
    !! | calc_x()            | calc_x()            | calc_x()               |
    !! | redistribute to sol | copy to sol         | redistribute to sol    |
    !! | interpolate to sol  |                     | interpolate to sol     |
    !!
    !! \return ierr
    integer function calc_norm_range(style,in_limits,eq_limits,X_limits,&
        &sol_limits,r_F_eq,r_F_X,r_F_sol,jq) result(ierr)
        
        character(*), parameter :: rout_name = 'calc_norm_range'
        
        ! input / output
        character(len=*), intent(in) :: style                                   !< style of calculation (PB3D: in, eq, X or sol; POST)
        integer, intent(inout), optional :: in_limits(2)                        !< min. and max. index of in grid
        integer, intent(inout), optional :: eq_limits(2)                        !< min. and max. index of eq grid for this process
        integer, intent(inout), optional :: X_limits(2)                         !< min. and max. index of X grid for this process
        integer, intent(inout), optional :: sol_limits(2)                       !< min. and max. index of sol grid for this process
        real(dp), intent(inout), optional :: r_F_eq(:)                          !< equilibrium r_F
        real(dp), intent(inout), allocatable, optional :: r_F_X(:)              !< perturbation r_F
        real(dp), intent(inout), optional :: r_F_sol(:)                         !< solution r_F
        real(dp), intent(in), optional :: jq(:)                                 !< q_saf (pol. flux) or rot_t (tor. flux) in total Flux variables
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! select depending on style
        select case (style)
            case ('PB3D_in')                                                    ! PB3D: in
                if (present(in_limits)) then
                    ierr = calc_norm_range_PB3D_in(in_limits)
                    CHCKERR('')
                else
                    ierr = 1
                    err_msg = 'for PB3D: in, in_limits need to be provided'
                    CHCKERR(err_msg)
                end if
            case ('PB3D_eq')                                                    ! PB3D: eq
                if (present(eq_limits)) then
                    call calc_norm_range_PB3D_eq(eq_limits)
                else
                    ierr = 1
                    err_msg = 'for PB3D: eq, eq_limits need to be provided'
                    CHCKERR(err_msg)
                end if
            case ('PB3D_X')                                                     ! PB3D: X
                if (present(eq_limits).and.present(X_limits).and.&
                    &present(r_F_eq).and.present(r_F_X).and.(present(jq))) then
                    ierr = calc_norm_range_PB3D_X(eq_limits,X_limits,&
                        &r_F_eq,r_F_X,jq)
                    CHCKERR('')
                else
                    ierr = 1
                    err_msg = 'for PB3D: X, eq_limits, X_limits, r_F_eq, &
                        &r_F_X and jq need to be provided'
                    CHCKERR(err_msg)
                end if
            case ('PB3D_sol')                                                   ! PB3D: sol
                if (present(sol_limits).and.present(r_F_sol)) then
                    ierr = calc_norm_range_PB3D_sol(sol_limits,r_F_sol)
                    CHCKERR('')
                else
                    ierr = 1
                    err_msg = 'for PB3D: sol, sol_limits and r_F_sol need to &
                        &be provided'
                    CHCKERR(err_msg)
                end if
            case ('POST')                                                       ! POST
                if (present(eq_limits) .and. present(X_limits) .and. &
                    &present(sol_limits) .and.  present(r_F_eq) .and. &
                    &present(r_F_X) .and.  present(r_F_sol)) then
                    call calc_norm_range_POST(eq_limits,X_limits,sol_limits,&
                        &r_F_eq,r_F_X,r_F_sol)
                else
                    ierr = 1
                    err_msg = 'for POST, eq_limits, X_limits, sol_limits, &
                        &r_F_eq, r_F_X and r_F_sol need to be provided'
                    CHCKERR(err_msg)
                end if
            case default
                ierr = 1
                err_msg = 'Incorrect style "'//trim(style)//'"'
                CHCKERR(err_msg)
        end select
    contains
        !> \public PB3D input version
        !!
        !! The normal range  is calculated by finding the tightest  range of the
        !! input variables encompassing the entire solution range. Additionally,
        !! the global max_flux variables are set as well..
        !!
        !! \note This is  the global, undivided range.  The division information
        !! is calculated  in 'eq_limits'. Furthermore, the  input variables have
        !! to be tabulated on the full  grid provided by the equilibrium code to
        !! be able to be used in PB3D and POST.
        integer function calc_norm_range_PB3D_in(in_limits) &
            &result(ierr)                                                       ! PB3D version for equilibrium grid
            use num_vars, only: use_pol_flux_E, use_pol_flux_F, eq_style, &
                &norm_disc_prec_eq
            use num_utilities, only: con2dis, dis2con, round_with_tol, spline
            use VMEC_vars, only: flux_p_V, flux_t_V
            use HELENA_vars, only: flux_p_H, flux_t_H
            use X_vars, only: min_r_sol, max_r_sol
            use eq_vars, only: max_flux_E, max_flux_F
            use grid_vars, only: n_r_in
            
            character(*), parameter :: rout_name = 'calc_norm_range_PB3D_in'
            
            ! input / output
            integer, intent(inout) :: in_limits(2)                              !< total min. and max. index of eq. grid for this process
            
            ! local variables
            real(dp), allocatable :: flux_F(:), flux_E(:)                       ! either pol. or tor. flux in F and E
            real(dp) :: tot_min_r_in_F_con                                      ! tot_min_r_in in continuous F coords.
            real(dp) :: tot_max_r_in_F_con                                      ! tot_max_r_in in continuous F coords.
            real(dp) :: tot_min_r_in_E_con(1)                                   ! tot_min_r_in in continuous E coords.
            real(dp) :: tot_max_r_in_E_con(1)                                   ! tot_max_r_in in continuous E coords.
            real(dp) :: tot_min_r_in_E_dis                                      ! tot_min_r_in in discrete E coords., unrounded
            real(dp) :: tot_max_r_in_E_dis                                      ! tot_max_r_in in discrete E coords., unrounded
            
            ! initialize ierr
            ierr = 0
            
            ! set up flux_F and flux_E,  depending on which equilibrium style is
            ! being used:
            !   1:  VMEC
            !   2:  HELENA
            allocate(flux_F(n_r_in),flux_E(n_r_in))
            select case (eq_style)
                case (1)                                                        ! VMEC
                    ! set up F flux
                    if (use_pol_flux_F) then
                        flux_F = flux_p_V(:,0)
                    else
                        flux_F = - flux_t_V(:,0)                                ! conversion VMEC LH -> RH coord. system
                    end if
                    ! set up E flux
                    if (use_pol_flux_E) then
                        flux_F = flux_p_V(:,0)
                    else
                        flux_E = flux_t_V(:,0)
                    end if
                case (2)                                                        ! HELENA
                    ! set up F flux
                    if (use_pol_flux_F) then
                        flux_F = flux_p_H(:,0)
                    else
                        flux_F = flux_t_H(:,0)
                    end if
                    ! set up E flux
                    if (use_pol_flux_E) then
                        flux_E = flux_p_H(:,0)
                    else
                        flux_E = flux_t_H(:,0)
                    end if
            end select
            
            ! set max_flux
            max_flux_E = flux_E(n_r_in)
            max_flux_F = flux_F(n_r_in)
            
            ! normalize flux_F and flux_E to (0..1) using max_flux
            flux_E = flux_E/max_flux_E
            flux_F = flux_F/max_flux_F
            
            ! get lower and upper bound of total solution range
            ! 1. start
            tot_min_r_in_F_con = min_r_sol
            tot_max_r_in_F_con = max_r_sol
            ! 2. round with tolerance
            ierr = round_with_tol(tot_min_r_in_F_con,0._dp,1._dp)
            CHCKERR('')
            ierr = round_with_tol(tot_max_r_in_F_con,0._dp,1._dp)
            CHCKERR('')
            ! 3. continuous E coords
            ierr = spline(flux_F,flux_E,[tot_min_r_in_F_con],&
                &tot_min_r_in_E_con,ord=min(3,norm_disc_prec_eq))
            ierr = spline(flux_F,flux_E,[tot_max_r_in_F_con],&
                &tot_max_r_in_E_con,ord=min(3,norm_disc_prec_eq))
            ! 4. round with tolerance
            ierr = round_with_tol(tot_min_r_in_E_con,minval(flux_E),&
                &maxval(flux_E))
            CHCKERR('')
            ierr = round_with_tol(tot_max_r_in_E_con,minval(flux_E),&
                &maxval(flux_E))
            CHCKERR('')
            ! 5. discrete E index, unrounded
            ierr = con2dis(tot_min_r_in_E_con(1),tot_min_r_in_E_dis,flux_E)
            CHCKERR('')
            ierr = con2dis(tot_max_r_in_E_con(1),tot_max_r_in_E_dis,flux_E)
            CHCKERR('')
            ! 6. discrete E index, rounded
            in_limits(1) = floor(tot_min_r_in_E_dis)
            in_limits(2) = ceiling(tot_max_r_in_E_dis)
        end function calc_norm_range_PB3D_in
        
        !> \public PB3D equilibrium version
        !!
        !! The normal range is calculated by dividing the full equilibrium range
        !! passed from the input phase between the processes.
        !!
        !! \note For  X_style 2 (solution) or  3 (optimized), at the  end of the
        !! equilibrium  phase, there  will  be  an exchange  of  data from  each
        !! process to each  process so that they all have  the tightest possible
        !! fit of data of the corresponding perturbation limits.
        subroutine calc_norm_range_PB3D_eq(eq_limits)                           ! PB3D version for equilibrium grid
            use num_vars, only: n_procs, rank, norm_disc_prec_eq
            use grid_vars, only: n_r_eq
            
            ! input / output
            integer, intent(inout) :: eq_limits(2)                              !< min. and max. index of eq. grid for this process
            
            ! set local equilibrium limits
            eq_limits(1) = nint(1 + 1._dp*rank*(n_r_eq-1)/n_procs)
            eq_limits(2) = nint(1 + (1._dp+rank)*(n_r_eq-1)/n_procs)
            
            ! increase lower limits for processes that are not first
            if (rank.gt.0) eq_limits(1) = eq_limits(1)+1
            
            ! ghost regions of width 2*norm_disc_prec_eq
            eq_limits(1) = max(eq_limits(1)-2*norm_disc_prec_eq,1)
            eq_limits(2) = min(eq_limits(2)+2*norm_disc_prec_eq,n_r_eq)
        end subroutine calc_norm_range_PB3D_eq
        
        !> \public PB3D perturbation version
        !!
        !! The normal range is determined according to X_grid style:
        !!  1. taken identical  to the equilibrium grid, which  means that after
        !!  the perturbation  phase the integrated tensorial  quantities have to
        !!  be  interpolated in  interp_v()  in driver_sol()  after having  been
        !!  redistributed at the start of the solution driver.
        !!  2.  taken identical  to  the  solution grid,  which  means that  the
        !!  the  equilibrium quantities  have  to be  interpolated in  calc_u(),
        !!  calc_kv() and calc_pv()  after having been redistributed  at the end
        !!  of the equilibrium driver.
        !!  3. by   considering  an   optimal  interpolation  extension  of  the
        !!  equilibrium range,  adding intermediairy points between  grid points
        !!  where  the  safety factor  changes  too  quickly, according  to  the
        !!  variable 'max_njq_change'. This uses \c prim_X.
        integer function calc_norm_range_PB3D_X(eq_limits,X_limits,r_F_eq,&
            &r_F_X,jq) result(ierr)                                             ! PB3D version for perturbation grid
            
            use num_vars, only: X_grid_style, max_njq_change, rank
            use grid_vars, only: n_r_eq, n_r_sol
            use X_vars, only: prim_X
            use grid_utilities, only: calc_eqd_grid
            use eq_vars, only: max_flux_F
            
            character(*), parameter :: rout_name = 'calc_norm_range_PB3D_X'
            
            ! input / output
            integer, intent(in) :: eq_limits(2)                                 !< min. and max. index of eq grid for this process
            integer, intent(inout) :: X_limits(2)                               !< min. and max. index of X grid for this process
            real(dp), intent(in) :: r_F_eq(:)                                   !< equilibrium r_F
            real(dp), intent(inout), allocatable :: r_F_X(:)                    !< perturbation r_F
            real(dp), intent(in) :: jq(:)                                       !< q_saf (pol. flux) or rot_t (tor. flux) in total Flux variables
            
            ! local variables
            integer :: kd                                                       ! counter
            integer, allocatable :: div(:)                                      ! number of extra divisions for this grid interval
            real(dp), allocatable :: r_F_plot(:,:)                              ! plot of r_F_eq and r_F_X
            
            ! initialize ierr
            ierr = 0
            
            select case (X_grid_style)
                case (1)                                                        ! equilibrium
                    ! copy equilibrium
                    allocate(r_F_X(n_r_eq))
                    r_F_X = r_F_eq
                    X_limits = eq_limits
                case (2)                                                        ! solution
                    ! calculate solution
                    allocate(r_F_X(n_r_sol))
                    ierr = calc_norm_range_PB3D_sol(sol_limits=X_limits,&
                        &r_F_sol=r_F_X)
                    CHCKERR('')
                case (3)                                                        ! enriched
                    ! decide extra divisions for each interval
                    allocate(div(size(r_F_eq)-1))
                    do kd = 1,size(r_F_eq)-1
                        div(kd) = floor(prim_X*(jq(kd+1)-jq(kd))/max_njq_change)
                    end do
                    
                    ! set up r_F_X with divisions
                    allocate(r_F_X(size(r_F_eq)+sum(div)))
                    do kd = 1,size(r_F_eq)-1
                        ierr = calc_eqd_grid(&
                            &r_F_X(kd+sum(div(1:kd-1)):kd+sum(div(1:kd))),&
                            &r_F_eq(kd),r_F_eq(kd+1),excl_last=.true.)
                        CHCKERR('')
                    end do
                    r_F_X(size(r_F_X)) = r_F_eq(size(r_F_eq))
                    
                    ! set up normal limits
                    X_limits(1) = eq_limits(1)+sum(div(1:eq_limits(1)-1))
                    X_limits(2) = eq_limits(2)+sum(div(1:eq_limits(2)-1))
                    
                    ! plot
                    if (rank.eq.0) then
                        allocate(r_F_plot(size(r_F_X),4))
                        r_F_plot(1:size(r_F_eq),1) = r_F_eq
                        r_F_plot(1:size(r_F_eq),3) = [(kd,kd=1,size(r_F_eq))]
                        r_F_plot(size(r_F_eq)+1:size(r_F_X),1) = &
                            &r_F_eq(size(r_F_eq))
                        r_F_plot(size(r_F_eq)+1:size(r_F_X),3) = size(r_F_eq)
                        r_F_plot(:,2) = r_F_X
                        r_F_plot(:,4) = [(kd,kd=1,size(r_F_X))]
                        r_F_plot(:,1:2) = r_F_plot(:,1:2)/max_flux_F*2*pi
                        call print_ex_2D(['eq','X '],'r_F_X',&
                            &r_F_plot(:,3:4),x=r_F_plot(:,1:2),draw=.false.)
                        call draw_ex(['eq','X '],'r_F_X',2,1,&
                            &.false.)
                    end if
            end select
        end function calc_norm_range_PB3D_X
        
        !> \public PB3D solution version
        !!
        !! The normal range is determined by simply dividing the solution range,
        !! a ghost range is required.
        !!
        !! By default, this routine uses  \c sol_n_procs processes, but this can
        !! be overruled.
        integer function calc_norm_range_PB3D_sol(sol_limits,r_F_sol,n_procs) &
            &result(ierr)                                                       ! PB3D version for solution grid
            use num_vars, only: sol_n_procs, rank, norm_disc_prec_sol
            use eq_vars, only: max_flux_F
            use X_vars, only: min_r_sol, max_r_sol
            use num_utilities, only: round_with_tol
            use grid_vars, only: n_r_sol
            
            character(*), parameter :: rout_name = 'calc_norm_range_PB3D_sol'
            
            ! input / output
            integer, intent(inout) :: sol_limits(2)                             !< min. and max. index of sol grid for this process
            real(dp), intent(inout) :: r_F_sol(:)                               !< solution r_F
            integer, intent(in), optional :: n_procs                            !< how many processes used
            
            ! local variables
            integer :: id                                                       ! counter
            integer, allocatable :: loc_n_r_sol(:)                              ! local nr. of normal points in solution grid
            integer :: n_procs_loc                                              ! local n_procs
            
            ! initialize ierr
            ierr = 0
            
            ! set local n_procs
            n_procs_loc = sol_n_procs
            if (present(n_procs)) n_procs_loc = n_procs
            
            ! set loc_n_r_sol
            allocate(loc_n_r_sol(n_procs_loc))
            loc_n_r_sol = n_r_sol/n_procs_loc                                   ! number of radial points on this processor
            loc_n_r_sol(1:mod(n_r_sol,n_procs_loc)) = &
                &loc_n_r_sol(1:mod(n_r_sol,n_procs_loc)) + 1                    ! add a mode to if there is a remainder
            
            ! set sol_limits
            if (rank.lt.n_procs_loc) then
                sol_limits = &
                    &[sum(loc_n_r_sol(1:rank))+1,sum(loc_n_r_sol(1:rank+1))]
                
                ! ghost regions of width 2*norm_disc_prec_sol
                sol_limits(1) = max(sol_limits(1)-norm_disc_prec_sol,1)
                sol_limits(2) = min(sol_limits(2)+norm_disc_prec_sol,n_r_sol)
            else
                sol_limits(1) = n_r_sol
                sol_limits(2) = n_r_sol
            end if
            
            ! calculate r_F_sol in range from min_r_sol to max_r_sol
            r_F_sol = [(min_r_sol + (id-1.0_dp)/(n_r_sol-1.0_dp)*&
                &(max_r_sol-min_r_sol),id=1,n_r_sol)]
            
            ! round with standard tolerance
            ierr = round_with_tol(r_F_sol,0.0_dp,1.0_dp)
            CHCKERR('')
            
            ! translate  to   the  real   normal  variable  in   range  from
            ! 0..flux/2pi
            r_F_sol = r_F_sol*max_flux_F/(2*pi)
        end function calc_norm_range_PB3D_sol
        
        !> \public POST version
        !!
        !! The normal range  is determined by simply dividing  a possible subset
        !! of the solution  range, indicated by \c min_r_plot  and \c max_r_sol,
        !! including  a  ghost range  and  getting  a bounding  equilibrium  and
        !! perturbation range.
        subroutine calc_norm_range_POST(eq_limits,X_limits,sol_limits,r_F_eq,&
            &r_F_X,r_F_sol)                                                     ! POST version
            
            use num_vars, only: n_procs, rank, norm_disc_prec_sol, &
                &min_r_plot, max_r_plot, X_grid_style
            use eq_vars, only: max_flux_F
            use grid_utilities, only: find_compr_range
            
            ! input / output
            integer, intent(inout) :: eq_limits(2)                              !< min. and max. index of eq grid for this process
            integer, intent(inout) :: X_limits(2)                               !< min. and max. index of X grid for this process
            integer, intent(inout) :: sol_limits(2)                             !< min. and max. index of sol grid for this process
            real(dp), intent(in) :: r_F_eq(:)                                   !< eq r_F
            real(dp), intent(in) :: r_F_X(:)                                    !< X r_F
            real(dp), intent(in) :: r_F_sol(:)                                  !< sol r_F
            
            ! local variables
            integer :: sol_limits_tot(2)                                        ! total solution limits
            integer :: n_r_sol                                                  ! total nr. of normal points in solution grid
            integer, allocatable :: loc_n_r_sol(:)                              ! local nr. of normal points in solution grid
            real(dp) :: min_sol, max_sol                                        ! min. and max. of r_F_sol in range of this process
            
            ! initialize ierr
            ierr = 0
            
            ! get min and max of solution range
            min_sol = max(minval(r_F_sol),min_r_plot*max_flux_F/(2*pi))
            max_sol = min(maxval(r_F_sol),max_r_plot*max_flux_F/(2*pi))
            
            ! find the solution index that comprises this range
            call find_compr_range(r_F_sol,[min_sol,max_sol],sol_limits_tot)
            
            ! initialize n_r_sol
            n_r_sol = sol_limits_tot(2)-sol_limits_tot(1)+1
            allocate(loc_n_r_sol(n_procs))
            
            ! divide the solution grid equally over all the processes
            loc_n_r_sol = n_r_sol/n_procs                                       ! number of radial points on this processor
            loc_n_r_sol(1:mod(n_r_sol,n_procs)) = &
                &loc_n_r_sol(1:mod(n_r_sol,n_procs)) + 1                        ! add a mode to if there is a remainder
            
            ! set sol_limits
            sol_limits = [sum(loc_n_r_sol(1:rank))+1,sum(loc_n_r_sol(1:rank+1))]
            sol_limits = sol_limits + sol_limits_tot(1) - 1
            if (rank.gt.0) sol_limits(1) = sol_limits(1)-norm_disc_prec_sol     ! ghost region for num. deriv.
            if (rank.lt.n_procs-1) sol_limits(2) = &
                &sol_limits(2)+norm_disc_prec_sol+1                             ! ghost region for num. deriv. and overlap for int_vol
            min_sol = minval(r_F_sol(sol_limits(1):sol_limits(2)))
            max_sol = maxval(r_F_sol(sol_limits(1):sol_limits(2)))
            
            ! determine eq_limits: smallest eq and X range comprising sol range
            call find_compr_range(r_F_eq,[min_sol,max_sol],eq_limits)
            select case (X_grid_style)
                case (1,3)                                                      ! equilibrium or enriched
                    call find_compr_range(r_F_X,[min_sol,max_sol],X_limits)
                case (2)                                                        ! solution
                    X_limits = sol_limits
            end select
        end subroutine calc_norm_range_POST
    end function calc_norm_range

    !> Sets up the equilibrium grid.
    !!
    !! The following  variables  are calculated:
    !!  - equilibrium variables (eq)
    !!  - perturbation variables (X)
    !!
    !! For the  implementation of the  equilibrium grid  the normal part  of the
    !! grid  is always  given by  the output  of the  equilibrium code,  but the
    !! angular part depends on the style:
    !!  -  VMEC:  The  output  is   analytic  in  the  angular  coordinates,  so
    !!  field-aligned coordinates are used.  Later, in calc_ang_grid_eq_b(), the
    !!  angular variables are calculated.
    !!  - HELENA: The output is given on a poloidal grid (no toroidal dependency
    !!  due to axisymmetry), which is used.
    !!
    !! \return ierr
    integer function setup_grid_eq(grid_eq,eq_limits) result(ierr)
        use num_vars, only: eq_style, eq_jobs_lims, eq_job_nr
        use grid_vars, only: n_r_eq, n_alpha
        use HELENA_vars, only: nchi, chi_H
        use rich_vars, only: n_par_X
        
        character(*), parameter :: rout_name = 'setup_grid_eq'
        
        ! input / output
        type(grid_type), intent(inout) :: grid_eq                               !< equilibrium grid
        integer, intent(in) :: eq_limits(2)                                     !< min. and max. index of eq grid of this process
        
        ! local variables
        integer :: id                                                           ! counter
        
        ! initialize ierr
        ierr = 0
        
        ! set up general equilibrium grid:
        ! choose which equilibrium style is being used:
        !   1:  VMEC
        !   2:  HELENA
        select case (eq_style)
            case (1)                                                            ! VMEC
                ! user output
                call writo('Field-aligned with '//trim(i2str(n_r_eq))//&
                    &' normal and '//trim(i2str(n_par_X))//' parallel points')
                
                ! create grid with eq_jobs_lims
                ierr = grid_eq%init([eq_jobs_lims(2,eq_job_nr)-&
                    &eq_jobs_lims(1,eq_job_nr)+1,n_alpha,n_r_eq],eq_limits)
                CHCKERR('')
            case (2)                                                            ! HELENA
                ! user output
                call writo('Identical to the equilibrium input grid')
                call writo(trim(i2str(n_r_eq))//' normal and '//&
                    &trim(i2str(nchi))//' angular points')
                call writo('Poloidal range '//trim(r2strt(chi_H(1)))//'..'//&
                    &trim(r2strt(chi_H(nchi))))
                call writo('Will be interpolated to a field-aligned grid &
                    &later')
                
                ! create grid
                ierr = grid_eq%init([nchi,1,n_r_eq],eq_limits)                  ! axisymmetric equilibrium
                CHCKERR('')
                
                ! copy angular grid from HELENA
                do id = 1,grid_eq%n(1)
                    grid_eq%theta_E(id,:,:) = chi_H(id)
                end do
                grid_eq%zeta_E = 0._dp
                
                ! convert to Flux coordinates (trivial)
                grid_eq%theta_F = grid_eq%theta_E
                grid_eq%zeta_F = grid_eq%zeta_E
        end select
    end function setup_grid_eq
    
    !> Sets up the field-aligned equilibrium grid.
    !!
    !! This serves  as a bridge  to the solution grid,  as it contains  the same
    !! normal coordinate  as the general  grid, but the angular  coordinates are
    !! defined by the solution grid.
    !!
    !! Optionally, only  half the  grid can  be calculated  (i.e. only  the even
    !! points), which is used for Richardson levels greater than 1.
    !!
    !! \note In contrast  to \c setup_grid_eq, the angular  coordinates are also
    !! calculated here.
    !!
    !! \return ierr
    integer function setup_grid_eq_B(grid_eq,grid_eq_B,eq,only_half_grid) &
        &result(ierr)
        use num_vars, only: eq_style, eq_jobs_lims, eq_job_nr
        use grid_vars, only: n_alpha
        
        character(*), parameter :: rout_name = 'setup_grid_eq_B'
        
        ! input / output
        type(grid_type), intent(inout) :: grid_eq                               !< general equilibrium grid
        type(grid_type), intent(inout) :: grid_eq_B                             !< field-aligned equilibrium grid
        type(eq_1_type), intent(in) :: eq                                       !< flux equilibrium variables
        logical, intent(in), optional :: only_half_grid                         !< calculate only half grid with even points
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! choose which equilibrium style is being used:
        !   1:  VMEC
        !   2:  HELENA
        select case (eq_style)
            case (1)                                                            ! VMEC
                ! the grid is already field aligned
                ierr = 1
                err_msg = 'The grid is already field-aligned for VMEC'
                CHCKERR(err_msg)
            case (2)                                                            ! HELENA
                ! create grid with eq_jobs_lims
                ierr = grid_eq_B%init([eq_jobs_lims(2,eq_job_nr)-&
                    &eq_jobs_lims(1,eq_job_nr)+1,n_alpha,grid_eq%n(3)],&
                    &[grid_eq%i_min,grid_eq%i_max])                             ! only one field line
                CHCKERR('')
                
                ! copy the normal coords.
                grid_eq_B%loc_r_E = grid_eq%loc_r_E
                grid_eq_B%loc_r_F = grid_eq%loc_r_F
                grid_eq_B%r_E = grid_eq%r_E
                grid_eq_B%r_F = grid_eq%r_F
                
                ! calculate the angular grid that follows the magnetic field
                ierr = calc_ang_grid_eq_B(grid_eq_B,eq,only_half_grid)
                CHCKERR('')
        end select
    end function setup_grid_eq_B
    
    !> Sets  up  the  general  perturbation  grid,  in  which  the  perturbation
    !! variables are calculated.
    !!
    !! For \c  X_grid_style 1, this grid  is identical to the  equilibrium grid,
    !! and  for \c  X_grid_style  2,it  has the  same  angular  extent but  with
    !! different normal points, indicated by the variable \c r_F_X.
    !!
    !! \return ierr
    integer function setup_grid_X(grid_eq,grid_X,r_F_X,X_limits) result(ierr)
        use num_vars, only: norm_disc_prec_X, X_grid_style
        use grid_vars, only: disc_type
        use grid_utilities, only: coord_F2E
        use num_utilities, only: spline
        
        character(*), parameter :: rout_name = 'setup_grid_X'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  !< equilibrium grid
        type(grid_type), intent(inout) :: grid_X                                !< perturbation grid
        real(dp), intent(in), allocatable :: r_F_X(:)                           !< points of perturbation grid
        integer, intent(in) :: X_limits(2)                                      !< min. and max. index of perturbation grid of this process
        
        ! local variables
        integer :: id, jd                                                       ! counters
        
        ! initialize ierr
        ierr = 0
        
        select case (X_grid_style)
            case (1)                                                            ! equilibrium
                ! X grid identical to equilibrium grid
                ierr = grid_eq%copy(grid_X)
                CHCKERR('')
            case (2,3)                                                          ! solution and enriched
                ! create grid
                ierr = grid_X%init([grid_eq%n(1:2),size(r_F_X)],X_limits)
                CHCKERR('')
                
                ! set Flux variables
                grid_X%r_F = r_F_X
                grid_X%loc_r_F = r_F_X(X_limits(1):X_limits(2))
                
                ! convert to Equilibrium variables
                ierr = coord_F2E(grid_eq,grid_X%r_F,grid_X%r_E,&
                    &r_F_array=grid_eq%r_F,r_E_array=grid_eq%r_E)
                CHCKERR('')
                ierr = coord_F2E(grid_eq,grid_X%loc_r_F,grid_X%loc_r_E,&
                    &r_F_array=grid_eq%r_F,r_E_array=grid_eq%r_E)
                CHCKERR('')
                
                ! interpolate
                do jd = 1,grid_eq%n(2)
                    do id = 1,grid_eq%n(1)
                        ierr = spline(grid_eq%loc_r_F,grid_eq%theta_E(id,jd,:),&
                            &grid_X%loc_r_F,grid_X%theta_E(id,jd,:),&
                            &ord=min(3,norm_disc_prec_X))
                        CHCKERR('')
                        ierr = spline(grid_eq%loc_r_F,grid_eq%zeta_E(id,jd,:),&
                            &grid_X%loc_r_F,grid_X%zeta_E(id,jd,:),&
                            &ord=min(3,norm_disc_prec_X))
                        CHCKERR('')
                        ierr = spline(grid_eq%loc_r_F,grid_eq%theta_F(id,jd,:),&
                            &grid_X%loc_r_F,grid_X%theta_F(id,jd,:),&
                            &ord=min(3,norm_disc_prec_X))
                        CHCKERR('')
                        ierr = spline(grid_eq%loc_r_F,grid_eq%zeta_F(id,jd,:),&
                            &grid_X%loc_r_F,grid_X%zeta_F(id,jd,:),&
                            &ord=min(3,norm_disc_prec_X))
                        CHCKERR('')
                    end do
                end do
        end select
    end function setup_grid_X
    
    !> Sets up  the general solution grid,  in which the solution  variables are
    !! calculated.
    !!
    !! For the  solution grid,  only one  parallel point  is used,  but possibly
    !! multiple geodesic points, equal to the number of field lines, \c n_alpha.
    !!
    !! \return ierr
    integer function setup_grid_sol(grid_eq,grid_X,grid_sol,r_F_sol,&
        &sol_limits) result(ierr)
        
        use num_vars, only: n_procs, X_grid_style
        use grid_utilities, only: coord_F2E
        
        character(*), parameter :: rout_name = 'setup_grid_sol'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  !< equilibrium grid
        type(grid_type), intent(in) :: grid_X                                   !< perturbation grid
        type(grid_type), intent(inout) :: grid_sol                              !< solution grid
        real(dp), intent(in) :: r_F_sol(:)                                      !< points of solution grid
        integer, intent(in) :: sol_limits(2)                                    !< min. and max. index of sol grid of this process
        
        ! initialize ierr
        ierr = 0
        
        select case (X_grid_style)
            case (1,3)                                                          ! equilibrium or enriched
                ! create grid
                ierr = grid_sol%init([0,0,size(r_F_sol)],sol_limits,&
                    &divided=n_procs.gt.1)
                CHCKERR('')
                
                ! set Flux variables
                grid_sol%r_F = r_F_sol
                grid_sol%loc_r_F = r_F_sol(sol_limits(1):sol_limits(2))
                
                ! convert to Equilibrium variables
                ierr = coord_F2E(grid_eq,grid_sol%r_F,grid_sol%r_E,&
                    &r_F_array=grid_eq%r_F,r_E_array=grid_eq%r_E)
                CHCKERR('')
                ierr = coord_F2E(grid_eq,grid_sol%loc_r_F,grid_sol%loc_r_E,&
                    &r_F_array=grid_eq%r_F,r_E_array=grid_eq%r_E)
                CHCKERR('')
            case (2)                                                            ! solution
                ! solution grid  identical to perturation  grid but with  only 1
                ! parallel point.
                ierr = grid_sol%init([0,0,grid_X%n(3)],&
                    &[grid_X%i_min,grid_X%i_max],grid_X%divided)
                CHCKERR('')
                grid_sol%r_F = grid_X%r_F
                grid_sol%loc_r_F = grid_X%r_F(sol_limits(1):sol_limits(2))
                grid_sol%r_E = grid_X%r_E
                grid_sol%loc_r_E = grid_X%r_E(sol_limits(1):sol_limits(2))
        end select
    end function setup_grid_sol
    
    !> Calculate equilibrium grid that follows magnetic field lines.
    !!
    !! This grid is  different from the equilibrium grid  from setup_grid_eq for
    !! HELENA, as the  latter is the output grid from  HELENA, which is situated
    !! in a single  poloidal cross-section, as opposed to  a really field-aligne
    !! grid.
    !!
    !! For VMEC, this is not used as the grid is field-aligned from the start.
    !! 
    !! \note
    !!  -# The end-points are included for  the grids in the parallel direction.
    !!  This is to facilitate working with the trapezoidal rule or Simpson's 3/8
    !!  rule for integration. This is \b NOT valid in general!
    !!  -# by  setting the flag \c  only_half_grid, only the even  points of the
    !!  parallel  grid are  calculated, which  is useful  for higher  Richardson
    !!  levels with VMEC so that only  new angular points are calculated and the
    !!  old ones reused.
    !!
    !! \return ierr
    integer function calc_ang_grid_eq_B(grid_eq,eq,only_half_grid) result(ierr)
        use num_vars, only: use_pol_flux_F, use_pol_flux_E, &
            &eq_style, tol_zero, eq_job_nr, eq_jobs_lims
        use grid_vars, only: min_par_X, max_par_X, alpha, n_alpha
        use eq_vars, only: max_flux_E
        use grid_utilities, only: coord_F2E, calc_eqd_grid, calc_n_par_X_rich
        use eq_utilities, only: print_info_eq
        use rich_vars, only: n_par_X
        use X_vars, only: min_r_sol, max_r_sol
#if ldebug
        use grid_utilities, only: coord_E2F
#endif
        
        character(*), parameter :: rout_name = 'calc_ang_grid_eq_B'
        
        ! input / output
        type(grid_type), intent(inout) :: grid_eq                               !< equilibrium grid of which to calculate angular part
        type(eq_1_type), intent(in), target :: eq                               !< flux equilibrium variables
        logical, intent(in), optional :: only_half_grid                         !< calculate only half grid with even points
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        real(dp), allocatable :: r_E_loc(:)                                     ! flux in Equilibrium coords.
        real(dp), pointer :: flux_F(:) => null()                                ! flux that the F uses as normal coord.
        real(dp), pointer :: flux_E(:) => null()                                ! flux that the E uses as normal coord.
        real(dp) :: r_F_factor, r_E_factor                                      ! mult. factors for r_F and r_E
        integer :: pmone                                                        ! plus or minus one
        integer :: id, jd, kd                                                   ! counters
        integer :: n_par_X_loc                                                  ! local n_par_X
        logical :: only_half_grid_loc                                           ! local only_half_grid
        real(dp), allocatable :: theta_F_loc(:,:,:)                             ! local theta_F
        real(dp), allocatable :: zeta_F_loc(:,:,:)                              ! local zeta_F
        
        ! initialize ierr
        ierr = 0
        
        call lvl_ud(1)
        
        ! set local only_half_grid
        only_half_grid_loc = .false.
        if (present(only_half_grid)) only_half_grid_loc = only_half_grid
        
        ! set local n_par_X_rich
        ierr = calc_n_par_X_rich(n_par_X_loc,only_half_grid_loc)
        CHCKERR('')
        
        ! user output
        call writo('for '//trim(i2str(grid_eq%n(3)))//&
            &' values on normal range '//trim(r2strt(min_r_sol))//'..'//&
            &trim(r2strt(max_r_sol)))
        call writo('for '//trim(i2str(n_par_X))//' values on parallel &
            &range '//trim(r2strt(min_par_X))//'pi..'//&
            &trim(r2strt(max_par_X))//'pi')
        call lvl_ud(1)
        if (only_half_grid_loc) call writo('for this Richardson level, only &
            &the even points are setup up')
        call print_info_eq(n_par_X_loc)
        call lvl_ud(-1)
        
        ! set up flux_E and plus minus one
        ! Note: this routine  is similar to calc_loc_r, but  that routine cannot
        ! be used here because it is possible that the FD quantities are not yet
        ! defined.
        ! choose which equilibrium style is being used:
        !   1:  VMEC
        !   2:  HELENA
        select case (eq_style)
            case (1)                                                            ! VMEC
                if (use_pol_flux_E) then                                        ! normalized poloidal flux
                    flux_E => eq%flux_p_E(:,0)
                else                                                            ! normalized toroidal flux
                    flux_E => eq%flux_t_E(:,0)
                end if
                r_E_factor = max_flux_E
                pmone = -1                                                      ! conversion VMEC LH -> RH coord. system
            case (2)                                                            ! HELENA
                flux_E => eq%flux_p_E(:,0)
                r_E_factor = 2*pi
                pmone = 1
        end select
        
        ! set up flux_F
        if (use_pol_flux_F) then                                                ! poloidal flux / 2pi
            flux_F => eq%flux_p_E(:,0)
            r_F_factor = 2*pi
        else                                                                    ! toroidal flux / 2pi
            flux_F => eq%flux_t_E(:,0)
            r_F_factor = pmone*2*pi                                             ! possible conversion VMEC LH -> RH coord. system
        end if
        
        ! set up parallel angle in  Flux coordinates on equidistant grid        
        ! and use this to calculate the other angle as well
        allocate(theta_F_loc(n_par_X,n_alpha,grid_eq%loc_n_r))
        allocate(zeta_F_loc(n_par_X,n_alpha,grid_eq%loc_n_r))
        if (use_pol_flux_F) then                                                ! parallel angle theta
            ierr = calc_eqd_grid(theta_F_loc,min_par_X*pi,&
                &max_par_X*pi,1,excl_last=.false.)                              ! first index corresponds to parallel angle
            CHCKERR('')
            do kd = 1,grid_eq%loc_n_r
                zeta_F_loc(:,:,kd) = pmone*eq%q_saf_E(kd,0)*&
                    &theta_F_loc(:,:,kd)
            end do
            do jd = 1,n_alpha
                zeta_F_loc(:,jd,:) = zeta_F_loc(:,jd,:) + alpha(jd)
            end do
        else                                                                    ! parallel angle zeta
            ierr = calc_eqd_grid(zeta_F_loc,min_par_X*pi,&
                &max_par_X*pi,1,excl_last=.false.)                              ! first index corresponds to parallel angle
            CHCKERR('')
            do kd = 1,grid_eq%loc_n_r
                theta_F_loc(:,:,kd) = pmone*eq%rot_t_E(kd,0)*&
                    &zeta_F_loc(:,:,kd)
            end do
            do jd = 1,n_alpha
                theta_F_loc(:,jd,:) = theta_F_loc(:,jd,:) - alpha(jd)
            end do
        end if
        
        ! set up grid_eq angular F coordinates
        if (only_half_grid_loc) then
            do id = eq_jobs_lims(1,eq_job_nr),eq_jobs_lims(2,eq_job_nr)
                grid_eq%theta_F(id-eq_jobs_lims(1,eq_job_nr)+1,:,:) = &
                    &theta_F_loc(2*id,:,:)
                grid_eq%zeta_F(id-eq_jobs_lims(1,eq_job_nr)+1,:,:) = &
                    &zeta_F_loc(2*id,:,:)
            end do
        else
            grid_eq%theta_F = theta_F_loc(eq_jobs_lims(1,eq_job_nr):&
                &eq_jobs_lims(2,eq_job_nr),:,:)
            grid_eq%zeta_F = zeta_F_loc(eq_jobs_lims(1,eq_job_nr):&
                &eq_jobs_lims(2,eq_job_nr),:,:)
        end if
        
        ! allocate local r_E
        allocate(r_E_loc(size(flux_F)))
        
        ! convert  Flux  coordinates  to  Equilibrium  coordinates  (use
        ! custom flux_E and  flux_F because the Flux  quantities are not
        ! yet calculated)
        call writo('convert F to E coordinates')
        call lvl_ud(1)
        ierr = coord_F2E(grid_eq,grid_eq%loc_r_F,grid_eq%theta_F,&
            &grid_eq%zeta_F,r_E_loc,grid_eq%theta_E,grid_eq%zeta_E,&
            &r_F_array=flux_F/r_F_factor,r_E_array=flux_E/r_E_factor)
        CHCKERR('')
        call lvl_ud(-1)
        
        ! test whether r_E_loc indeed corresponds to loc_r_E of eq grid
        if (maxval(abs(grid_eq%loc_r_E-r_E_loc)).gt.10*tol_zero) then
            ierr = 1
            err_msg = 'loc_r_E of equilibrium grid is not recovered'
            CHCKERR(err_msg)
        end if
        
#if ldebug
        if (debug_calc_ang_grid_eq_B) then
            ! test whether F variables recovered
            deallocate(theta_F_loc,zeta_F_loc)
            allocate(theta_F_loc(grid_eq%n(1),grid_eq%n(2),grid_eq%loc_n_r))
            allocate(zeta_F_loc(grid_eq%n(1),grid_eq%n(2),grid_eq%loc_n_r))
            ierr = coord_E2F(grid_eq,grid_eq%loc_r_E,grid_eq%theta_E,&
                &grid_eq%zeta_E,grid_eq%loc_r_F,theta_F_loc,zeta_F_loc,&
                &r_E_array=flux_E/r_E_factor,r_F_array=flux_F/r_F_factor)
            CHCKERR('')
            
            ! plot difference
            call plot_diff_HDF5(grid_eq%theta_F,theta_F_loc,'TEST_theta_F',&
                &grid_eq%n,[0,0,grid_eq%i_min-1],'test whether F variable are &
                &recovered',output_message=.true.)
            call plot_diff_HDF5(grid_eq%zeta_F,zeta_F_loc,'TEST_zeta_F',&
                &grid_eq%n,[0,0,grid_eq%i_min-1],'test whether F variable are &
                &recovered',output_message=.true.)
        end if
#endif
        
        ! deallocate local variables
        deallocate(r_E_loc)
        nullify(flux_F,flux_E)
        
        call lvl_ud(-1)
    end function calc_ang_grid_eq_B
    
    !> Redistribute a grid to match the normal distribution of solution grid.
    !! 
    !! The routine first calculates the smallest eq range that comprises the sol
    !! range.  Then, it  gets the  lowest equilibrium  limits able  to setup  an
    !! output grid that starts at index 1. After determining the output grid, it
    !! then sends the variables to their new processes using MPI.
    !!
    !! \note
    !!  -# Only the  Flux variables are saved.
    !!  -# the redistributed grid has trimmed  outer limits, i.e. it starts at 1
    !!  and ends at the upper limit of  the last process. This can be turned off
    !!  optionally using \c no_outer_trim.
    !!
    !! \return ierr
    integer function redistribute_output_grid(grid,grid_out,no_outer_trim) &
        &result(ierr)
        use grid_utilities, only: find_compr_range
        use MPI_utilities, only: redistribute_var, get_ser_var
        use grid_vars, only: n_r_sol
        use num_vars, only: n_procs
        
        character(*), parameter :: rout_name = 'redistribute_output_grid'
        
        ! input / output
        type(grid_type), intent(in) :: grid                                     !< equilibrium grid variables
        type(grid_type), intent(inout) :: grid_out                              !< redistributed equilibrium grid variables
        logical, intent(in), optional :: no_outer_trim                          !< do not trim the outer limits
        
        ! local variables
        integer :: id                                                           ! counter
        integer :: eq_limits(2)                                                 ! normal limits for equilibrium variables
        integer :: sol_limits(2)                                                ! normal limits for perturbation variables
        integer :: n_out(3)                                                     ! n of grid_out
        integer :: i_lim_tot(2)                                                 ! total limits of grid
        integer :: i_lim_out(2)                                                 ! limits of grid_out
        integer :: lims(2), lims_dis(2)                                         ! limits and distributed limits, taking into account the angular extent
        real(dp), allocatable :: r_F_sol(:)                                     ! perturbation r_F
        real(dp), allocatable :: temp_var(:)                                    ! temporary variable
        integer, allocatable :: temp_lim(:)                                     ! temporary limit
        integer, allocatable :: eq_limits_tot(:,:)                              ! total equilibrium limits
        logical :: no_outer_trim_loc                                            ! local no_outer_trim
        
        ! initialize ierr
        ierr = 0
        
        ! calculate normal range for solution and save in perturbation variables
        allocate(r_F_sol(n_r_sol))
        ierr = calc_norm_range('PB3D_sol',sol_limits=sol_limits,r_F_sol=r_F_sol)
        CHCKERR('')
        
        ! determine eq_limits: smallest eq range comprising X range
        call find_compr_range(grid%r_F,r_F_sol(sol_limits),eq_limits)
        
        ! get lowest equilibrium limits to be  able to setup an output grid that
        ! starts at index 1
        allocate(eq_limits_tot(2,n_procs))
        do id = 1,2
            ierr = get_ser_var(eq_limits(id:id),temp_lim,scatter=.true.)
            CHCKERR('')
            eq_limits_tot(id,:) = temp_lim
        end do
        
        ! set up redistributed grid
        no_outer_trim_loc = .false.
        if (present(no_outer_trim)) no_outer_trim_loc = no_outer_trim
        if (no_outer_trim_loc) then
            n_out = grid%n
            i_lim_tot = [1,n_out(3)]
            i_lim_out = eq_limits
        else
            n_out(1:2) = grid%n(1:2)
            n_out(3) = eq_limits_tot(2,n_procs)-eq_limits_tot(1,1)+1
            i_lim_tot = [eq_limits_tot(1,1),eq_limits_tot(2,n_procs)]
            i_lim_out = eq_limits-eq_limits_tot(1,1)+1
        end if
        ierr = grid_out%init(n_out,i_lim_out)
        CHCKERR('')
        
        ! set up limits taking into account angular extent and temporary var
        lims(1) = product(grid%n(1:2))*(grid%i_min-1)+1
        lims(2) = product(grid%n(1:2))*grid%i_max
        lims_dis(1) = product(grid%n(1:2))*(eq_limits(1)-1)+1
        lims_dis(2) = product(grid%n(1:2))*eq_limits(2)
        allocate(temp_var(lims_dis(2)-lims_dis(1)+1))
        
        ! copy total variables
        ! r_F
        grid_out%r_F = grid%r_F(i_lim_tot(1):i_lim_tot(2))
        ! r_E
        grid_out%r_E = grid%r_E(i_lim_tot(1):i_lim_tot(2))
        
        ! distribute local variables
        ! theta_F
        ierr = redistribute_var(reshape(grid%theta_F,[size(grid%theta_F)]),&
            &temp_var,lims,lims_dis)
        CHCKERR('')
        grid_out%theta_F = reshape(temp_var,shape(grid_out%theta_F))
        ! zeta_F
        ierr = redistribute_var(reshape(grid%zeta_F,[size(grid%zeta_F)]),&
            &temp_var,lims,lims_dis)
        CHCKERR('')
        grid_out%zeta_F = reshape(temp_var,shape(grid_out%zeta_F))
        ! theta_E
        ierr = redistribute_var(reshape(grid%theta_E,[size(grid%theta_E)]),&
            &temp_var,lims,lims_dis)
        CHCKERR('')
        grid_out%theta_E = reshape(temp_var,shape(grid_out%theta_E))
        ! zeta_E
        ierr = redistribute_var(reshape(grid%zeta_E,[size(grid%zeta_E)]),&
            &temp_var,lims,lims_dis)
        CHCKERR('')
        grid_out%zeta_E = reshape(temp_var,shape(grid_out%zeta_E))
        ! loc_r_F and loc_R_E
        if (grid_out%divided) then
            grid_out%loc_r_F = grid_out%r_F(grid_out%i_min:grid_out%i_max)
            grid_out%loc_r_E = grid_out%r_E(grid_out%i_min:grid_out%i_max)
        end if
    end function redistribute_output_grid
    
    !> Plots the grid in real 3-D space.
    !!
    !! This creates an animation that can be used by ParaView or VisIt.
    !!
    !! The equilibrium grid should contain the fieldline-oriented angles with \c
    !! ang_1 the parallel angle and \c ang_2 the field line label.
    !!
    !! \see See grid_vars.grid_type for a discussion on \c ang_1 and \c ang_2.
    !!
    !! \note
    !!  -# This procedure  does not use \c n_theta_plot and  \c n_zeta_plot from
    !!  num_vars, but instead temporarily overwrites them with its own, since it
    !!  is suposed to be 3-D also in the axisymmetric case.
    !!  -# The implementation is currently very slow.
    !!
    !! \return ierr
    integer function magn_grid_plot(grid) result(ierr)
        use num_vars, only: rank, no_plots, n_theta_plot, n_zeta_plot, &
            &eq_style, min_theta_plot, max_theta_plot, min_zeta_plot, &
            &max_zeta_plot
        use grid_utilities, only: trim_grid, extend_grid_F, calc_XYZ_grid
        use VMEC_utilities, only: calc_trigon_factors
        
        character(*), parameter :: rout_name = 'magn_grid_plot'
        
        ! input / output
        type(grid_type), intent(in) :: grid                                     !< fieldline-oriented equilibrium grid
        
        ! local variables
        real(dp), allocatable :: X_1(:,:,:), Y_1(:,:,:), Z_1(:,:,:)             ! X, Y and Z of surface in Axisymmetric coordinates
        real(dp), allocatable :: X_2(:,:,:), Y_2(:,:,:), Z_2(:,:,:)             ! X, Y and Z of magnetic field lines in Axisymmetric coordinates
        real(dp), pointer :: X_1_tot(:,:,:) => null()                           ! total X
        real(dp), pointer :: Y_1_tot(:,:,:) => null()                           ! total Y
        real(dp), pointer :: Z_1_tot(:,:,:) => null()                           ! total Z
        real(dp), pointer :: X_2_tot(:,:,:) => null()                           ! total X
        real(dp), pointer :: Y_2_tot(:,:,:) => null()                           ! total Y
        real(dp), pointer :: Z_2_tot(:,:,:) => null()                           ! total Z
        type(grid_type) :: grid_ext                                             ! angularly extended grid
        type(grid_type) :: grid_plot                                            ! grid for plotting
        integer :: n_theta_plot_old                                             ! backup of n_theta_plot
        integer :: n_zeta_plot_old                                              ! backup of n_zeta_plot
        real(dp) :: min_theta_plot_old, max_theta_plot_old                      ! backup of min and max_theta_plot
        real(dp) :: min_zeta_plot_old, max_zeta_plot_old                        ! backup of min and max_zeta_plot
        character(len=max_str_ln) :: anim_name                                  ! name of animation
        
        ! initialize ierr
        ierr = 0
        
        ! bypass plots if no_plots
        if (no_plots) return
        
        call writo('Plotting magnetic field and flux surfaces')
        
        call lvl_ud(1)
        
        ! 0. set up variables
        
        ! save n_theta_plot and n_zeta_plot and change them
        n_theta_plot_old = n_theta_plot
        n_zeta_plot_old = n_zeta_plot
        min_theta_plot_old = min_theta_plot
        max_theta_plot_old = max_theta_plot
        min_zeta_plot_old = min_zeta_plot
        max_zeta_plot_old = max_zeta_plot
        n_theta_plot = 40
        n_zeta_plot = 160
        min_theta_plot = 1
        max_theta_plot = 3
        min_zeta_plot = 0
        max_zeta_plot = 2
        
        ! extend grid
        ierr = extend_grid_F(grid,grid_ext,grid_eq=grid)
        CHCKERR('')
        
        ! restore n_theta_plot and n_zeta_plot
        n_theta_plot = n_theta_plot_old
        n_zeta_plot = n_zeta_plot_old
        min_theta_plot = min_theta_plot_old
        max_theta_plot = max_theta_plot_old
        min_zeta_plot = min_zeta_plot_old
        max_zeta_plot = max_zeta_plot_old
        
        ! trim extended grid into plot grid
        ierr = trim_grid(grid_ext,grid_plot)
        CHCKERR('')
        
        ! if VMEC, calculate trigonometric factors of plot grid
        if (eq_style.eq.1) then
            ierr = calc_trigon_factors(grid_plot%theta_E,grid_plot%zeta_E,&
                &grid_plot%trigon_factors)
            CHCKERR('')
        end if
        
        ! set animation name
        anim_name = 'Magnetic field in flux surfaces'
        
        ! 1. plot flux surfaces
        call writo('writing flux surfaces')
        
        ! calculate X_1,Y_1 and Z_1
        allocate(X_1(grid_plot%n(1),grid_plot%n(2),grid_plot%loc_n_r))
        allocate(Y_1(grid_plot%n(1),grid_plot%n(2),grid_plot%loc_n_r))
        allocate(Z_1(grid_plot%n(1),grid_plot%n(2),grid_plot%loc_n_r))
        ierr = calc_XYZ_grid(grid,grid_plot,X_1,Y_1,Z_1)
        CHCKERR('')
        
        ! dealloc grid
        call grid_plot%dealloc()
        
        ! get pointers to full X, Y and Z
        ! The reason for this is that the plot  is not as simple as usual, so no
        ! divided  plots  are used,  and  also  efficiency  is not  the  biggest
        ! priority. Therefore,  all the plotting  of the  local is handled  by a
        ! single process, the master.
        call get_full_XYZ(X_1,Y_1,Z_1,X_1_tot,Y_1_tot,Z_1_tot,'flux surfaces')
        
        ! 2. plot field lines
        call writo('writing field lines')
        
        ! trim grid into plot grid
        ierr = trim_grid(grid,grid_plot)
        CHCKERR('')
        
        ! if VMEC, calculate trigonometric factors of plot grid
        if (eq_style.eq.1) then
            ierr = calc_trigon_factors(grid_plot%theta_E,grid_plot%zeta_E,&
                &grid_plot%trigon_factors)
            CHCKERR('')
        end if
        
        ! calculate X_2,Y_2 and Z_2
        allocate(X_2(grid_plot%n(1),grid_plot%n(2),grid_plot%loc_n_r))
        allocate(Y_2(grid_plot%n(1),grid_plot%n(2),grid_plot%loc_n_r))
        allocate(Z_2(grid_plot%n(1),grid_plot%n(2),grid_plot%loc_n_r))
        ierr = calc_XYZ_grid(grid,grid_plot,X_2,Y_2,Z_2)
        CHCKERR('')
        
        ! dealloc grids
        call grid_plot%dealloc()
        call grid_ext%dealloc()
        
        ! get pointers to full X, Y and Z
        ! The reason for this is that the plot  is not as simple as usual, so no
        ! divided  plots  are used,  and  also  efficiency  is not  the  biggest
        ! priority. Therefore,  all the plotting  of the  local is handled  by a
        ! single process, the master.
        call get_full_XYZ(X_2,Y_2,Z_2,X_2_tot,Y_2_tot,Z_2_tot,'field lines')
        
        ierr = magn_grid_plot_HDF5(X_1_tot,X_2_tot,Y_1_tot,Y_2_tot,&
            &Z_1_tot,Z_2_tot,anim_name)
        CHCKERR('')
        
        ! clean up
        deallocate(X_1,Y_1,Z_1)
        deallocate(X_2,Y_2,Z_2)
        nullify(X_1_tot,Y_1_tot,Z_1_tot)
        nullify(X_2_tot,Y_2_tot,Z_2_tot)
        
        call lvl_ud(-1)
        
        call writo('Done plotting magnetic field and flux surfaces')
    contains
        ! get pointer to full plotting variables X, Y and Z
        !> \private
        subroutine get_full_XYZ(X,Y,Z,X_tot,Y_tot,Z_tot,merge_name)
            use MPI_utilities, only: get_ser_var
            
            ! input / output
            real(dp), intent(in), target :: X(:,:,:), Y(:,:,:), Z(:,:,:)        ! X, Y and Z of either flux surfaces or magnetic field lines
            real(dp), intent(inout), pointer :: X_tot(:,:,:)                    ! pointer to full X
            real(dp), intent(inout), pointer :: Y_tot(:,:,:)                    ! pointer to full Y
            real(dp), intent(inout), pointer :: Z_tot(:,:,:)                    ! pointer to full Z
            character(len=*) :: merge_name                                      ! name of variable to be merged
            
            ! local variables
            real(dp), allocatable :: ser_XYZ_loc(:)                             ! serial copy of XYZ_loc
            integer, allocatable :: tot_dim(:)                                  ! total dimensions for plot 
            
            ! merge plots for flux surfaces if more than one process
            if (grid%divided) then                                              ! merge local plots
                ! user output
                call writo('Merging local plots for '//merge_name)
                
                ! get total dimension
                ierr = get_ser_var([size(X,3)],tot_dim)
                CHCKERR('')
                
                ! allocate total X, Y and Z (should have same sizes)
                if (rank.eq.0) then
                    allocate(X_tot(size(X,1),size(X,2),sum(tot_dim)))
                    allocate(Y_tot(size(X,1),size(X,2),sum(tot_dim)))
                    allocate(Z_tot(size(X,1),size(X,2),sum(tot_dim)))
                end if
                
                allocate(ser_XYZ_loc(size(X(:,:,1))*sum(tot_dim)))
                ierr = get_ser_var(reshape(X,[size(X)]),ser_XYZ_loc)
                CHCKERR('')
                if (rank.eq.0) X_tot = reshape(ser_XYZ_loc,shape(X_tot))
                ierr = get_ser_var(reshape(Y,[size(Y)]),ser_XYZ_loc)
                CHCKERR('')
                if (rank.eq.0) Y_tot = reshape(ser_XYZ_loc,shape(Y_tot))
                ierr = get_ser_var(reshape(Z,[size(Z)]),ser_XYZ_loc)
                CHCKERR('')
                if (rank.eq.0) Z_tot = reshape(ser_XYZ_loc,shape(Z_tot))
            else                                                                ! just point
                X_tot => X
                Y_tot => Y
                Z_tot => Z
            end if
        end subroutine get_full_XYZ
        
        ! Plot with HDF5
        !> \private
        integer function magn_grid_plot_HDF5(X_1,X_2,Y_1,Y_2,Z_1,Z_2,&
            &anim_name) result(ierr)
            use HDF5_ops, only: open_HDF5_file, add_HDF5_item, &
                &print_HDF5_top, print_HDF5_geom, print_HDF5_3D_data_item, &
                &print_HDF5_grid, close_HDF5_file
            use HDF5_vars, only: dealloc_XML_str, &
                & XML_str_type, HDF5_file_type
            use rich_vars, only: rich_lvl
            use grid_vars, only: n_alpha
            
            character(*), parameter :: rout_name = 'magn_grid_plot_HDF5'
          
            ! input / output
            real(dp), intent(in), pointer :: X_1(:,:,:), Y_1(:,:,:), Z_1(:,:,:) ! X, Y and Z of surface in Axisymmetric coordinates
            real(dp), intent(in), pointer :: X_2(:,:,:), Y_2(:,:,:), Z_2(:,:,:) ! X, Y and Z of magnetic field lines in Axisymmetric coordinates
            character(len=*), intent(in) :: anim_name                           ! name of animation
            
            ! local variables
            character(len=max_str_ln) :: plot_title(2)                          ! plot title for flux surface and for field lines
            character(len=max_str_ln) :: file_name                              ! name of file
            integer :: id, jd                                                   ! counter
            type(HDF5_file_type) :: file_info                                   ! file info
            type(XML_str_type), allocatable :: grids(:)                         ! grid in spatial collection (flux surface, field line)
            type(XML_str_type), allocatable :: space_col_grids(:)               ! grid with space collection at different times
            type(XML_str_type) :: time_col_grid                                 ! grid with time collection
            type(XML_str_type) :: top                                           ! topology
            type(XML_str_type) :: XYZ(3)                                        ! data items for geometry
            type(XML_str_type) :: geom                                          ! geometry
            integer :: loc_dim(3,2)                                             ! local dimensions (flux surfaces, field lines)
            integer :: n_r                                                      ! nr. of normal points
            
            ! initialize ierr
            ierr = 0
            
            ! only master
            if (rank.eq.0) then
                ! user output
                call writo('Drawing animation with HDF5')
                call lvl_ud(1)
                
                ! set up loc_dim and n_r
                loc_dim(:,1) = [size(X_1,1),size(X_1,2),1]
                loc_dim(:,2) = [size(X_2,1),1,1]
                n_r = size(X_1,3) - 1                                           ! should be same for all other X_i, Y_i and Z_i
                
                ! set up plot titles and file name
                plot_title(1) = 'Magnetic Flux Surfaces'
                plot_title(2) = 'Magnetic Field Lines'
                file_name = 'field_lines_in_flux_surfaces_R_'//&
                    &trim(i2str(rich_lvl))
                
                ! open HDF5 file
                ierr = open_HDF5_file(file_info,trim(file_name),&
                    &descr=anim_name,ind_plot=.true.)
                CHCKERR('')
                
                ! create grid  for flux surface and all field  lines and one for
                ! time collection
                allocate(grids(1+n_alpha))
                allocate(space_col_grids(n_r))
                
                ! loop over all normal points
                do id = 1,n_r
                    ! A. start with flux surface
                    
                    ! print topology
                    call print_HDF5_top(top,2,loc_dim(:,1))
                    
                    ! print data item for X
                    ierr = print_HDF5_3D_data_item(XYZ(1),file_info,&
                        &'X_surf_'//trim(i2str(id)),X_1(:,:,id:id),loc_dim(:,1))
                    CHCKERR('')
                    
                    ! print data item for Y
                    ierr = print_HDF5_3D_data_item(XYZ(2),file_info,&
                        &'Y_surf_'//trim(i2str(id)),Y_1(:,:,id:id),loc_dim(:,1))
                    CHCKERR('')
                    
                    ! print data item for Z
                    ierr = print_HDF5_3D_data_item(XYZ(3),file_info,&
                        &'Z_surf_'//trim(i2str(id)),Z_1(:,:,id:id),loc_dim(:,1))
                    CHCKERR('')
                    
                    ! print geometry with X, Y and Z data item
                    call print_HDF5_geom(geom,2,XYZ,.true.)
                    
                    ! create a grid with the topology and the geometry
                    ierr = print_HDF5_grid(grids(1),plot_title(1),1,&
                        &grid_top=top,grid_geom=geom,reset=.true.)
                    CHCKERR('')
                    
                    ! B. end with magnetic field
                    
                    do jd = 1,n_alpha                                           ! iterate over all field lines
                        ! print topology
                        call print_HDF5_top(top,2,loc_dim(:,2))
                        
                        ! print data item for X of this field line
                        ierr = print_HDF5_3D_data_item(XYZ(1),file_info,&
                            &'X_field_'//trim(i2str(id))//'_'//trim(i2str(jd)),&
                            &X_2(:,jd:jd,id+1:id+1),loc_dim(:,2))
                        CHCKERR('')
                        
                        ! print data item for Y of this field line
                        ierr = print_HDF5_3D_data_item(XYZ(2),file_info,&
                            &'Y_field_'//trim(i2str(id))//'_'//trim(i2str(jd)),&
                            &Y_2(:,jd:jd,id+1:id+1),loc_dim(:,2))
                        CHCKERR('')
                        
                        ! print data item for Z of this field line
                        ierr = print_HDF5_3D_data_item(XYZ(3),file_info,&
                            &'Z_field_'//trim(i2str(id))//'_'//trim(i2str(jd)),&
                            &Z_2(:,jd:jd,id+1:id+1),loc_dim(:,2))
                        CHCKERR('')
                        
                        ! print geometry with X, Y and Z data item
                        call print_HDF5_geom(geom,2,XYZ,.true.)
                        
                        ! create a grid with the topology and the geometry
                        ierr = print_HDF5_grid(grids(jd+1),plot_title(2),1,&
                            &grid_top=top,grid_geom=geom,reset=.true.)
                        CHCKERR('')
                    end do
                    
                    ! C.  merge flux  surface and magnetic  field in  to spatial
                    ! collection
                    ierr = print_HDF5_grid(space_col_grids(id),&
                        &'spatial collection',3,grid_time=id*1._dp,&
                        &grid_grids=grids,reset=.true.)
                    CHCKERR('')
                end do
                
                ! create grid collection from individual grids and reset them
                ierr = print_HDF5_grid(time_col_grid,'time collection',2,&
                    &grid_grids=space_col_grids,reset=.true.)
                CHCKERR('')
                
                ! add collection grid to HDF5 file and reset it
                ierr = add_HDF5_item(file_info,time_col_grid,reset=.true.)
                CHCKERR('')
                
                ! close HDF5 file
                ierr = close_HDF5_file(file_info)
                CHCKERR('')
                
                call lvl_ud(-1)
                
                ! clean up
                do id = 1,2
                    call dealloc_XML_str(grids(id))
                end do
                call dealloc_XML_str(space_col_grids)
                call dealloc_XML_str(time_col_grid)
                call dealloc_XML_str(top)
                do id = 1,3
                    call dealloc_XML_str(XYZ(id))
                end do
                call dealloc_XML_str(geom)
            end if
        end function magn_grid_plot_HDF5
    end function magn_grid_plot
    
    !> Print grid variables to an output file.
    !!
    !! If \c  rich_lvl is  provided, <tt>_R_[rich_lvl]</tt>  is appended  to the
    !! data name if it is > 0.
    !!
    !! Optionally, it  can be specified  that this  is a divided  parallel grid,
    !! corresponding to the variable \c eq_jobs_lims with index \c eq_job_nr. In
    !! this case,  the total grid  size is adjusted to  the one specified  by \c
    !! eq_jobs_lims and the grid is written as a subset.
    !!
    !! \note <tt>grid_</tt> is added in front the data_name.
    !!
    !! \return ierr
    integer function print_output_grid(grid,grid_name,data_name,rich_lvl,&
        &par_div,remove_previous_arrs) result(ierr)
        use num_vars, only: PB3D_name, eq_jobs_lims, eq_job_nr
        use HDF5_ops, only: print_HDF5_arrs
        use HDF5_vars, only: var_1D_type, &
            &max_dim_var_1D
        use grid_utilities, only: trim_grid
        
        character(*), parameter :: rout_name = 'print_output_grid'
        
        ! input / output
        type(grid_type), intent(in) :: grid                                     !< grid variables
        character(len=*), intent(in) :: grid_name                               !< name to display
        character(len=*), intent(in) :: data_name                               !< name under which to store
        integer, intent(in), optional :: rich_lvl                               !< Richardson level to reconstruct
        logical, intent(in), optional :: par_div                                !< is a parallely divided grid
        logical, intent(in), optional :: remove_previous_arrs                   !< remove previous variables if present
        
        ! local variables
        integer :: n_tot(3)                                                     ! total n
        integer :: par_id(2)                                                    ! local parallel interval
        type(grid_type) :: grid_trim                                            ! trimmed grid
        type(var_1D_type), allocatable, target :: grid_1D(:)                    ! 1D equivalent of grid variables
        type(var_1D_type), pointer :: grid_1D_loc => null()                     ! local element in grid_1D
        integer :: id                                                           ! counter
        logical :: par_div_loc                                                  ! local par_div
        
        ! initialize ierr
        ierr = 0
        
        ! user output
        call writo('Write '//trim(grid_name)//' grid variables to output &
            &file')
        call lvl_ud(1)
        
        ! trim grid
        ierr = trim_grid(grid,grid_trim)
        CHCKERR('')
        
        ! set local par_div
        par_div_loc = .false.
        if (present(par_div)) par_div_loc = par_div
        
        ! set total n and parallel interval
        n_tot = grid_trim%n
        par_id = [1,n_tot(1)]
        if (grid_trim%n(1).gt.0 .and. par_div_loc) then                         ! total grid includes all equilibrium jobs
            n_tot(1) = maxval(eq_jobs_lims)-minval(eq_jobs_lims)+1
            par_id = eq_jobs_lims(:,eq_job_nr)
        end if
        
        ! Set up the 1D equivalents of the perturbation variables
        allocate(grid_1D(max_dim_var_1D))
        
        ! set up variables grid_1D
        id = 1
        
        ! n
        grid_1D_loc => grid_1D(id); id = id+1
        grid_1D_loc%var_name = 'n'
        allocate(grid_1D_loc%tot_i_min(1),grid_1D_loc%tot_i_max(1))
        allocate(grid_1D_loc%loc_i_min(1),grid_1D_loc%loc_i_max(1))
        grid_1D_loc%tot_i_min = [1]
        grid_1D_loc%tot_i_max = [3]
        grid_1D_loc%loc_i_min = [1]
        grid_1D_loc%loc_i_max = [3]
        allocate(grid_1D_loc%p(3))
        grid_1D_loc%p = n_tot
        
        ! r_F
        grid_1D_loc => grid_1D(id); id = id+1
        grid_1D_loc%var_name = 'r_F'
        allocate(grid_1D_loc%tot_i_min(1),grid_1D_loc%tot_i_max(1))
        allocate(grid_1D_loc%loc_i_min(1),grid_1D_loc%loc_i_max(1))
        grid_1D_loc%tot_i_min = [1]
        grid_1D_loc%tot_i_max = [n_tot(3)]
        grid_1D_loc%loc_i_min = [grid_trim%i_min]
        grid_1D_loc%loc_i_max = [grid_trim%i_max]
        allocate(grid_1D_loc%p(size(grid_trim%loc_r_F)))
        grid_1D_loc%p = grid_trim%loc_r_F
        
        ! r_E
        grid_1D_loc => grid_1D(id); id = id+1
        grid_1D_loc%var_name = 'r_E'
        allocate(grid_1D_loc%tot_i_min(1),grid_1D_loc%tot_i_max(1))
        allocate(grid_1D_loc%loc_i_min(1),grid_1D_loc%loc_i_max(1))
        grid_1D_loc%tot_i_min = [1]
        grid_1D_loc%tot_i_max = [n_tot(3)]
        grid_1D_loc%loc_i_min = [grid_trim%i_min]
        grid_1D_loc%loc_i_max = [grid_trim%i_max]
        allocate(grid_1D_loc%p(size(grid_trim%loc_r_E)))
        grid_1D_loc%p = grid_trim%loc_r_E
        
        ! only for 3D grids
        if (product(grid%n(1:2)).ne.0) then
            ! theta_F
            grid_1D_loc => grid_1D(id); id = id+1
            grid_1D_loc%var_name = 'theta_F'
            allocate(grid_1D_loc%tot_i_min(3),grid_1D_loc%tot_i_max(3))
            allocate(grid_1D_loc%loc_i_min(3),grid_1D_loc%loc_i_max(3))
            grid_1D_loc%tot_i_min = [1,1,1]
            grid_1D_loc%tot_i_max = n_tot
            grid_1D_loc%loc_i_min = [par_id(1),1,grid_trim%i_min]
            grid_1D_loc%loc_i_max = [par_id(2),n_tot(2),grid_trim%i_max]
            allocate(grid_1D_loc%p(size(grid_trim%theta_F)))
            grid_1D_loc%p = reshape(grid_trim%theta_F,[size(grid_trim%theta_F)])
            
            ! theta_E
            grid_1D_loc => grid_1D(id); id = id+1
            grid_1D_loc%var_name = 'theta_E'
            allocate(grid_1D_loc%tot_i_min(3),grid_1D_loc%tot_i_max(3))
            allocate(grid_1D_loc%loc_i_min(3),grid_1D_loc%loc_i_max(3))
            grid_1D_loc%tot_i_min = [1,1,1]
            grid_1D_loc%tot_i_max = n_tot
            grid_1D_loc%loc_i_min = [par_id(1),1,grid_trim%i_min]
            grid_1D_loc%loc_i_max = [par_id(2),n_tot(2),grid_trim%i_max]
            allocate(grid_1D_loc%p(size(grid_trim%theta_E)))
            grid_1D_loc%p = reshape(grid_trim%theta_E,[size(grid_trim%theta_E)])
            
            ! zeta_F
            grid_1D_loc => grid_1D(id); id = id+1
            grid_1D_loc%var_name = 'zeta_F'
            allocate(grid_1D_loc%tot_i_min(3),grid_1D_loc%tot_i_max(3))
            allocate(grid_1D_loc%loc_i_min(3),grid_1D_loc%loc_i_max(3))
            grid_1D_loc%tot_i_min = [1,1,1]
            grid_1D_loc%tot_i_max = n_tot
            grid_1D_loc%loc_i_min = [par_id(1),1,grid_trim%i_min]
            grid_1D_loc%loc_i_max = [par_id(2),n_tot(2),grid_trim%i_max]
            allocate(grid_1D_loc%p(size(grid_trim%zeta_F)))
            grid_1D_loc%p = reshape(grid_trim%zeta_F,[size(grid_trim%zeta_F)])
            
            ! zeta_E
            grid_1D_loc => grid_1D(id); id = id+1
            grid_1D_loc%var_name = 'zeta_E'
            allocate(grid_1D_loc%tot_i_min(3),grid_1D_loc%tot_i_max(3))
            allocate(grid_1D_loc%loc_i_min(3),grid_1D_loc%loc_i_max(3))
            grid_1D_loc%tot_i_min = [1,1,1]
            grid_1D_loc%tot_i_max = n_tot
            grid_1D_loc%loc_i_min = [par_id(1),1,grid_trim%i_min]
            grid_1D_loc%loc_i_max = [par_id(2),n_tot(2),grid_trim%i_max]
            allocate(grid_1D_loc%p(size(grid_trim%zeta_E)))
            grid_1D_loc%p = reshape(grid_trim%zeta_E,[size(grid_trim%zeta_E)])
        end if
        
        ! write
        ierr = print_HDF5_arrs(grid_1D(1:id-1),PB3D_name,&
            &'grid_'//trim(data_name),rich_lvl=rich_lvl,&
            &ind_print=.not.grid_trim%divided,&
            &remove_previous_arrs=remove_previous_arrs)
        CHCKERR('')
        
        ! clean up
        call grid_trim%dealloc()
        
        ! clean up
        nullify(grid_1D_loc)
        
        ! user output
        call lvl_ud(-1)
    end function print_output_grid
end module grid_ops
