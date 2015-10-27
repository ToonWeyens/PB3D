!------------------------------------------------------------------------------!
!   Operations that have to do with the grids and different coordinate systems !
!------------------------------------------------------------------------------!
module grid_ops
#include <PB3D_macros.h>
    use str_ops
    use output_ops
    use messages
    use num_vars, only: dp, pi, max_str_ln
    use grid_vars, only: grid_type
    use eq_vars, only: eq_type
    use met_vars, only: met_type

    implicit none
    private
    public coord_F2E, coord_E2F, calc_XYZ_grid, calc_norm_range, &
        &calc_eqd_grid, calc_ang_grid_eq, calc_ang_grid_eq_B, plot_grid_real, &
        &extend_grid_E, setup_grid_eq, setup_and_calc_grid_eq_B, &
        &setup_and_calc_grid_sol, get_norm_interp_data, calc_int_vol, &
        &trim_grid, untrim_grid
#if ldebug
    public debug_calc_ang_grid_eq_B, debug_get_norm_interp_data, &
        &debug_calc_int_vol
#endif
    
    interface coord_F2E
        module procedure coord_F2E_r, coord_F2E_rtz
    end interface
    interface coord_E2F
        module procedure coord_E2F_r, coord_E2F_rtz
    end interface
    interface calc_eqd_grid
        module procedure calc_eqd_grid_1D, calc_eqd_grid_3D
    end interface
    
    ! global variables
#if ldebug
    logical :: debug_calc_ang_grid_eq_B = .false.                               ! plot debug information for calc_ang_grid_eq_B
    logical :: debug_get_norm_interp_data = .false.                             ! plot debug information for get_norm_interp_data
    logical :: debug_calc_int_vol = .false.                                     ! plot debug information for calc_int_vol
#endif
    
contains
    ! Calculates normal range for the equilibrium grid and/or the solution grid.
    ! For PB3D:
    !   - if  'eq_limits'  is  provided,  they  are  calculated by  finding  the
    !   tightest  equilibrium grid  points  that encompass  the entire  solution
    !   range. 'r_F_eq' is ignored.
    !   - if 'sol_limits'  is provided, the solution range is  divided under the
    !   processes and  'r_F_eq' is filled. Note  that 'r_F_eq' has to  be of the
    !   correct size.
    ! For POST:
    !   -  The  solution  range  is calculated  and  the  tightest  encompassing
    !   equilibrium range is found. Also 'r_F_eq' and 'r_F_sol' are provided and
    !   necessary.
    integer function calc_norm_range(eq_limits,sol_limits,r_F_eq,r_F_sol) &
        &result(ierr)
        use num_vars, only: prog_style
        
        character(*), parameter :: rout_name = 'calc_norm_range'
        
        ! input / output
        real(dp), intent(inout), optional :: r_F_eq(:)                          ! equilibrium r_F
        real(dp), intent(inout), optional :: r_F_sol(:)                         ! solution r_F
        integer, intent(inout), optional :: eq_limits(2)                        ! min. and max. index of eq grid for this process
        integer, intent(inout), optional :: sol_limits(2)                       ! min. and max. index of sol grid for this process
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! select depending on program style
        select case (prog_style)
            case(1)                                                             ! PB3D
                if (present(eq_limits)) then
                    ierr = calc_norm_range_PB3D_eq(eq_limits)
                    CHCKERR('')
                end if
                if (present(sol_limits).and.present(r_F_sol)) then
                    ierr = calc_norm_range_PB3D_sol(sol_limits,r_F_sol)
                    CHCKERR('')
                end if
            case(2)                                                             ! PB3D post-processing
                if (.not.present(eq_limits) .or. .not.present(sol_limits) .or. &
                    &.not.present(r_F_eq) .or. .not. present(r_F_sol)) then
                    ierr = 1
                    CHCKERR('Incorrect variables provided.')
                end if
                ierr = calc_norm_range_POST(eq_limits,sol_limits,r_F_eq,r_F_sol)
                CHCKERR('')
            case default
                err_msg = 'No program style associated with '//&
                    &trim(i2str(prog_style))
                ierr = 1
                CHCKERR(err_msg)
        end select
    contains
        ! The normal range  is calculated by finding  tightest equilibrium range
        ! encompassing  the  entire  solution  range,  simply  dividing  it  and
        ! including a ghost region.
        integer function calc_norm_range_PB3D_eq(eq_limits) result(ierr)         ! PB3D version for equilibrium grid
            use num_vars, only: n_procs, rank, use_pol_flux_E, &
                &use_pol_flux_F, eq_style, tol_norm_r, norm_disc_prec_eq
            use utilities, only: con2dis, dis2con, calc_int, interp_fun, &
                &calc_deriv, round_with_tol
            use VMEC, only: flux_t_V, Dflux_t_V, rot_t_V
            use HELENA, only: flux_p_H, qs
            use X_vars, only: min_r_sol, max_r_sol
            use grid_vars, only: n_r_eq
            
            character(*), parameter :: rout_name = 'calc_norm_range_PB3D_eq'
            
            ! input / output
            integer, intent(inout) :: eq_limits(2)                              ! min. and max. index of eq. grid for this process
            
            ! local variables
            real(dp), allocatable :: flux_F(:), flux_E(:)                       ! either pol. or tor. flux in F and E
            real(dp), allocatable :: Dflux_p_H(:)                               ! normal derivative of flux_p_H
            real(dp) :: tot_min_r_eq_F_con                                      ! tot_min_r_eq in continuous F coords.
            real(dp) :: tot_max_r_eq_F_con                                      ! tot_max_r_eq in continuous F coords.
            real(dp) :: tot_min_r_eq_E_con                                      ! tot_min_r_eq in continuous E coords.
            real(dp) :: tot_max_r_eq_E_con                                      ! tot_max_r_eq in continuous E coords.
            real(dp) :: tot_min_r_eq_E_dis                                      ! tot_min_r_eq in discrete E coords., unrounded
            real(dp) :: tot_max_r_eq_E_dis                                      ! tot_max_r_eq in discrete E coords., unrounded
            integer :: tot_eq_limits(2)                                         ! total min. and max. index of eq. grid
            character(len=max_str_ln) :: err_msg                                ! error message
            
            ! initialize ierr
            ierr = 0
            
            ! set up flux_F and flux_E,  depending on which equilibrium style is
            ! being used:
            !   1:  VMEC
            !   2:  HELENA
            allocate(flux_F(n_r_eq),flux_E(n_r_eq))
            select case (eq_style)
                case (1)                                                        ! VMEC
                    ! set up F flux
                    if (use_pol_flux_F) then
                        ierr = calc_int(-rot_t_V*Dflux_t_V,&
                            &1.0_dp/(n_r_eq-1.0_dp),flux_F)                     ! conversion VMEC LH -> RH coord. system
                        CHCKERR('')
                    else
                        flux_F = flux_t_V
                    end if
                    ! set up E flux
                    if (use_pol_flux_E) then
                    write(*,*) '!!! CHANGED rot_t_V to -rot_t_V BUT NOT SURE !!'
                        ierr = calc_int(-rot_t_V*Dflux_t_V,&
                            &1.0_dp/(n_r_eq-1.0_dp),flux_E)                     ! conversion VMEC LH -> RH coord. system
                        CHCKERR('')
                    else
                        flux_E = flux_t_V
                    end if
                case (2)                                                        ! HELENA
                    ! calculate normal derivative of flux_p_H
                    allocate(Dflux_p_H(n_r_eq))
                    ierr = calc_deriv(flux_p_H,Dflux_p_H,flux_p_H,1,&
                        &norm_disc_prec_eq)
                    CHCKERR('')
                    ! set up F flux
                    if (use_pol_flux_F) then
                        flux_F = flux_p_H
                    else
                        ierr = calc_int(qs*Dflux_p_H,flux_p_H,flux_F)
                        CHCKERR('')
                    end if
                    ! set up E flux
                    if (use_pol_flux_E) then
                        flux_E = flux_p_H
                    else
                        ierr = calc_int(qs*Dflux_p_H,flux_p_H,flux_E)
                        CHCKERR('')
                    end if
                case default
                    err_msg = 'No equilibrium style associated with '//&
                        &trim(i2str(eq_style))
                    ierr = 1
                    CHCKERR(err_msg)
            end select
            
            ! normalize flux_F and flux_E to (0..1)
            ! Note: max_flux  is not yet  determined, so use  flux_F(n_r_eq) and
            ! flux_E(n_r_eq), knowing that this is the full global range
            flux_F = flux_F/flux_F(n_r_eq)
            flux_E = flux_E/flux_E(n_r_eq)
            
            ! get lower and upper bound of total solution range
            ! 1. include tolerance
            tot_min_r_eq_F_con = max(min_r_sol-tol_norm_r,0._dp)
            tot_max_r_eq_F_con = min(max_r_sol+tol_norm_r,1._dp)
            ! 2. round with tolerance
            ierr = round_with_tol(tot_min_r_eq_F_con,0._dp,1._dp)
            CHCKERR('')
            ierr = round_with_tol(tot_max_r_eq_F_con,0._dp,1._dp)
            CHCKERR('')
            ! 3. continuous E coords
            ierr = interp_fun(tot_min_r_eq_E_con,flux_E,tot_min_r_eq_F_con,&
                &flux_F)
            CHCKERR('')
            ierr = interp_fun(tot_max_r_eq_E_con,flux_E,tot_max_r_eq_F_con,&
                &flux_F)
            CHCKERR('')
            ! 4. discrete E index, unrounded
            ierr = con2dis(tot_min_r_eq_E_con,tot_min_r_eq_E_dis,flux_E)
            CHCKERR('')
            ierr = con2dis(tot_max_r_eq_E_con,tot_max_r_eq_E_dis,flux_E)
            CHCKERR('')
            ! 5. discrete E index, rounded
            tot_eq_limits(1) = floor(tot_min_r_eq_E_dis)
            tot_eq_limits(2) = ceiling(tot_max_r_eq_E_dis)
            
            ! set local equilibrium limits
            eq_limits(1) = nint(tot_eq_limits(1) + 1._dp*rank*&
                &(tot_eq_limits(2)-tot_eq_limits(1))/n_procs)
            eq_limits(2) = nint(tot_eq_limits(1) + (1._dp+rank)*&
                &(tot_eq_limits(2)-tot_eq_limits(1))/n_procs)
            
            ! increase lower limits for processes that are not first
            if (rank.gt.0) eq_limits(1) = eq_limits(1)+1
            
            ! ghost regions of width 2*norm_disc_prec_eq
            eq_limits(1) = &
                &max(eq_limits(1)-2*norm_disc_prec_eq,tot_eq_limits(1))
            eq_limits(2) = &
                &min(eq_limits(2)+2*norm_disc_prec_eq,tot_eq_limits(2))
            
            ! clean up
            deallocate(flux_F,flux_E)
        end function calc_norm_range_PB3D_eq
        
        ! The normal range  is determined by simply dividing  the solution range
        ! and including a ghost region.
        integer function calc_norm_range_PB3D_sol(sol_limits,r_F_sol) &
            &result(ierr)                                                       ! PB3D version for solution grid
            use num_vars, only: n_procs, rank, norm_disc_prec_sol, &
                &use_pol_flux_F
            use eq_vars, only: max_flux_p_F, max_flux_t_F
            use X_vars, only: min_r_sol, max_r_sol
            use utilities, only: round_with_tol
            
            character(*), parameter :: rout_name = 'calc_norm_range_PB3D_sol'
            
            ! input / output
            integer, intent(inout) :: sol_limits(2)                             ! min. and max. index of sol grid for this process
            real(dp), intent(inout) :: r_F_sol(:)                               ! solution r_F
            
            ! local variables
            integer :: n_r_sol                                                  ! total nr. of normal points in solution grid
            integer :: id                                                       ! counter
            integer, allocatable :: loc_n_r_sol(:)                              ! local nr. of normal points in solution grid
            real(dp) :: max_flux_F                                              ! either max_flux_p_F or max_flux_t_F
            
            ! initialize ierr
            ierr = 0
            
            ! set loc_n_r_sol
            n_r_sol = size(r_F_sol)
            allocate(loc_n_r_sol(n_procs))
            loc_n_r_sol = n_r_sol/n_procs                                       ! number of radial points on this processor
            loc_n_r_sol(1:mod(n_r_sol,n_procs)) = &
                &loc_n_r_sol(1:mod(n_r_sol,n_procs)) + 1                        ! add a mode to if there is a remainder
            
            ! set sol_limits
            sol_limits = [sum(loc_n_r_sol(1:rank))+1,sum(loc_n_r_sol(1:rank+1))]
            if (rank.gt.0) sol_limits(1) = sol_limits(1)-norm_disc_prec_sol     ! ghost region for num. deriv.
            if (rank.lt.n_procs-1) sol_limits(2) = &
                &sol_limits(2)+norm_disc_prec_sol+1                             ! ghost region for num. deriv. and overlap for int_vol
            
            
            ! set up max_flux
            if (use_pol_flux_F) then
                max_flux_F = max_flux_p_F
            else
                max_flux_F = max_flux_t_F
            end if
            
            ! calculate r_F_sol in range from 0 to 1
            r_F_sol = [(min_r_sol + (id-1.0_dp)/(n_r_sol-1.0_dp)*&
                &(max_r_sol-min_r_sol),id=1,n_r_sol)]
            
            ! round with standard tolerance
            ierr = round_with_tol(r_F_sol,0.0_dp,1.0_dp)
            CHCKERR('')
            
            ! translate to the real normal variable in range from 0..flux/2pi
            r_F_sol = r_F_sol*max_flux_F/(2*pi)
        end function calc_norm_range_PB3D_sol
        
        ! The normal range is determined  by simply dividing the solution range,
        ! including a ghost range and getting a bounding equilibrium range.
        integer function calc_norm_range_POST(eq_limits,sol_limits,r_F_eq,&
            &r_F_sol) result(ierr)                                              ! POST version
            use num_vars, only: n_procs, rank, norm_disc_prec_sol
            use utilities, only: con2dis, dis2con, calc_int, interp_fun, &
                &calc_deriv, round_with_tol
            
            character(*), parameter :: rout_name = 'calc_norm_range_POST'
            
            ! input / output
            integer, intent(inout) :: eq_limits(2), sol_limits(2)               ! min. and max. index of eq and sol grid for this process
            real(dp), intent(in) :: r_F_eq(:), r_F_sol(:)                       ! equilibrium and solution r_F
            
            ! local variables
            integer :: n_r_eq, n_r_sol                                          ! total nr. of normal points in eq and solution grid
            integer, allocatable :: loc_n_r_eq(:), loc_n_r_sol(:)               ! local nr. of normal points in eq and solution grid
            integer :: id                                                       ! counter
            real(dp) :: min_sol, max_sol                                        ! min. and max. of r_F_sol in range of this process
            real(dp), parameter :: tol = 1.E-6                                  ! tolerance for grids
            character(len=max_str_ln) :: err_msg                                ! error message
            
            ! initialize ierr
            ierr = 0
            
            ! initialize n_r_eq and n_r_sol
            n_r_eq = size(r_F_eq)
            n_r_sol = size(r_F_sol)
            allocate(loc_n_r_eq(n_procs),loc_n_r_sol(n_procs))
            
            ! divide the solution grid equally over all the processes
            loc_n_r_sol = n_r_sol/n_procs                                       ! number of radial points on this processor
            loc_n_r_sol(1:mod(n_r_sol,n_procs)) = &
                &loc_n_r_sol(1:mod(n_r_sol,n_procs)) + 1                        ! add a mode to if there is a remainder
            
            ! set sol_limits
            sol_limits = [sum(loc_n_r_sol(1:rank))+1,sum(loc_n_r_sol(1:rank+1))]
            if (rank.gt.0) sol_limits(1) = sol_limits(1)-norm_disc_prec_sol     ! ghost region for num. deriv.
            if (rank.lt.n_procs-1) sol_limits(2) = &
                &sol_limits(2)+norm_disc_prec_sol+1                             ! ghost region for num. deriv. and overlap for int_vol
            min_sol = minval(r_F_sol(sol_limits(1):sol_limits(2)))
            max_sol = maxval(r_F_sol(sol_limits(1):sol_limits(2)))
            
            ! determine eq_limits: smallest eq range comprising sol range
            eq_limits = [0,size(r_F_eq)+1]
            if (r_F_eq(1).lt.r_F_eq(size(r_F_eq))) then                         ! ascending r_F_eq
                do id = 1,size(r_F_eq)
                    if (r_F_eq(id).le.min_sol+tol) eq_limits(1) = id            ! move lower limit up
                    if (r_F_eq(size(r_F_eq)-id+1).ge.max_sol-tol) &
                        &eq_limits(2) = size(r_F_eq)-id+1                       ! move upper limit down
                end do
            else                                                                ! descending r_F_eq
                do id = 1,size(r_F_eq)
                    if (r_F_eq(size(r_F_eq)-id+1).le.min_sol+tol) &
                        &eq_limits(1) = size(r_F_eq)-id+1                       ! move lower limit up
                    if (r_F_eq(id).ge.max_sol-tol) eq_limits(2) = id            ! move upper limit down
                end do
            end if
            
            ! check if valid limits found
            if (eq_limits(1).lt.1 .or. eq_limits(2).gt.size(r_F_eq)) then
                ierr = 1
                err_msg = 'Solution range not contained in equilibrium range'
                CHCKERR(err_msg)
            end if
        end function calc_norm_range_POST
    end function calc_norm_range

    ! Sets up the general equilibrium grid, in which the following variables are
    ! calculated:
    !   - equilibrium variables (eq)
    !   - metric variables (met)
    !   - perturbation variables (X)
    ! For the implementation of the equilibrium grid the normal part of the grid
    ! is always  given by the  output of the  equilibrium code, but  the angular
    ! part depends on the style:
    !   -  VMEC:  The  output  is   analytic  in  the  angular  coordinates,  so
    !   field-oriented coordinates are used.
    !   - HELENA: The output is given on a poloidal grid (no toroidal dependency
    !   due to axisymmetry), so
    ! Here,  the  normal   variables  are  set  up  in  this   grid.  Later,  in
    ! calc_ang_grid_eq, the angular variables are calculated.
    ! Note that the general equilibrium grid using HELENA does not depend on the
    ! perturbation, and  therefore could theoretically be  calculated only once,
    ! before entering the Richardson loop.  However, to keep the structure clear
    ! this is not done and the grid is recalculated for every step.
    integer function setup_grid_eq(grid_eq,eq_limits) result(ierr)
        use num_vars, only: eq_style
        use grid_vars, only: create_grid, &
            &n_par_X, n_r_eq
        use HELENA, only: nchi
        
        character(*), parameter :: rout_name = 'setup_grid_eq'
        
        ! input / output
        type(grid_type), intent(inout) :: grid_eq                               ! equilibrium grid
        integer, intent(in) :: eq_limits(2)                                     ! min. and max. index of eq grid of this process
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! set up general equilibrium grid:
        ! choose which equilibrium style is being used:
        !   1:  VMEC
        !   2:  HELENA
        select case (eq_style)
            case (1)                                                            ! VMEC
                ierr = create_grid(grid_eq,[n_par_X,1,n_r_eq],eq_limits)        ! only one field line
                CHCKERR('')
            case (2)                                                            ! HELENA
                ierr = create_grid(grid_eq,[nchi,1,n_r_eq],eq_limits)           ! axisymmetric equilibrium
                CHCKERR('')
            case default
                err_msg = 'No equilibrium style associated with '//&
                    &trim(i2str(eq_style))
                ierr = 1
                CHCKERR(err_msg)
        end select
    end function setup_grid_eq
    
    ! Sets up  the field-aligned equilibrium grid,  which serves as a  bridge to
    ! the  solution grid,  as  it contains  the same  normal  coordinate as  the
    ! general  grid, but  the angular  coordinates are  defined by  the solution
    ! grid.
    ! In contrast to setup_grid_eq,  the angular coordinates are also calculated
    ! here.
    integer function setup_and_calc_grid_eq_B(grid_eq,grid_eq_B,eq) &
        &result(ierr)
        use num_vars, only: eq_style, plot_grid
        use grid_vars, only: create_grid, &
            &n_par_X
        
        character(*), parameter :: rout_name = 'setup_and_calc_grid_eq_B'
        
        ! input / output
        type(grid_type), intent(inout), target :: grid_eq                       ! general equilibrium grid
        type(grid_type), intent(inout), pointer :: grid_eq_B                    ! field-aligned equilibrium grid
        type(eq_type), intent(in) :: eq                                         ! equilibrium containing the angular grid
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! choose which equilibrium style is being used:
        !   1:  VMEC
        !   2:  HELENA
        select case (eq_style)
            case (1)                                                            ! VMEC
                ! user output
                call writo('Equilibrium grid is already field-aligned')
                call lvl_ud(1)
                
                ! the general equilibrium grid is already field-aligned
                grid_eq_B => grid_eq
            case (2)                                                            ! HELENA
                ! user output
                call writo('Start determining field-aligned grid in '//&
                    &trim(i2str(n_par_X))//' parallel and '//&
                    &trim(i2str(grid_eq%n(3)))//' normal points')
                call lvl_ud(1)
                
                ! create the grid
                allocate(grid_eq_B)
                ierr = create_grid(grid_eq_B,[n_par_X,1,grid_eq%n(3)],&
                    &[grid_eq%i_min,grid_eq%i_max])
                CHCKERR('')
                
                ! copy the normal coords.
                grid_eq_B%loc_r_E = grid_eq%loc_r_E
                grid_eq_B%loc_r_F = grid_eq%loc_r_F
                grid_eq_B%r_E = grid_eq%r_E
                grid_eq_B%r_F = grid_eq%r_F
                
                ! calculate the angular grid that follows the magnetic field
                ierr = calc_ang_grid_eq_B(grid_eq_B,eq)
                CHCKERR('')
            case default
                err_msg = 'No equilibrium style associated with '//&
                    &trim(i2str(eq_style))
                ierr = 1
                CHCKERR(err_msg)
        end select
        
        ! plot grid if requested
        if (plot_grid) then
            ierr = plot_grid_real(grid_eq_B)
            CHCKERR('')
        else
            call writo('Magnetic grid plot not requested')
        end if
        
        ! user output
        call lvl_ud(-1)
        call writo('Field-aligned equilibrium grid set up')
    end function setup_and_calc_grid_eq_B
    
    ! Sets  up the general  solution grid, in  which the solution  variables are
    ! calculated.
    integer function setup_and_calc_grid_sol(grid_eq,grid_sol,eq,r_F_sol,&
        &sol_limits) result(ierr)
        use grid_vars, only: create_grid
        
        character(*), parameter :: rout_name = 'setup_and_calc_grid_sol'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid
        type(grid_type), intent(inout) :: grid_sol                              ! solution grid
        type(eq_type), intent(in) :: eq                                         ! equilibrium variables
        real(dp), intent(in) :: r_F_sol(:)                                      ! points of solution grid
        integer, intent(in) :: sol_limits(2)                                    ! min. and max. index of sol grid of this process
        
        ! initialize ierr
        ierr = 0
        
        ! create grid
        ierr = create_grid(grid_sol,size(r_F_sol),sol_limits)
        CHCKERR('')
        
        ! set Flux variables
        grid_sol%r_F = r_F_sol
        grid_sol%loc_r_F = r_F_sol(sol_limits(1):sol_limits(2))
        
        ! convert to Equilibrium variables
        ierr = coord_F2E(grid_eq,eq,grid_sol%r_F,grid_sol%r_E,&
            &r_F_array=grid_eq%r_F,r_E_array=grid_eq%r_E)
        CHCKERR('')
        ierr = coord_F2E(grid_eq,eq,grid_sol%loc_r_F,grid_sol%loc_r_E,&
            &r_F_array=grid_eq%r_F,r_E_array=grid_eq%r_E)
        CHCKERR('')
    end function setup_and_calc_grid_sol
    
    ! Calculate the angular equilibrium grid.
    ! As explained in setup_grid_eq, the angular part of the equilibrium grid is
    ! determined by the equilibrium style.
    ! The variable  use_pol_flux determines  whether theta (.true.)  or zeta
    ! (.false.) is used as the parallel variable.
    integer function calc_ang_grid_eq(grid_eq,eq) result(ierr)
        use num_vars, only: eq_style
        use HELENA, only: chi_H
        
        character(*), parameter :: rout_name = 'calc_ang_grid_eq'
        
        ! input / output
        type(grid_type), intent(inout) :: grid_eq                               ! general equilibrium grid
        type(eq_type), intent(in) :: eq                                         ! equilibrium containing the angular grid
        
        ! local variables
        integer :: jd, kd                                                       ! counters
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! choose which equilibrium style is being used:
        !   1:  VMEC
        !   2:  HELENA
        select case (eq_style)
            case (1)                                                            ! VMEC
                ! calculate the angular grid that follows the magnetic field
                ierr = calc_ang_grid_eq_B(grid_eq,eq)
                CHCKERR('')
            case (2)                                                            ! HELENA
                ! calculate the angular grid from HELENA
                do kd = 1,grid_eq%loc_n_r
                    do jd = 1,grid_eq%n(2)
                        grid_eq%theta_E(:,jd,kd) = chi_H
                    end do
                end do
                grid_eq%zeta_E = 0._dp
                
                ! convert to Flux coordinates (trivial)
                grid_eq%theta_F = grid_eq%theta_E
                grid_eq%zeta_F = grid_eq%zeta_E
            case default
                err_msg = 'No equilibrium style associated with '//&
                    &trim(i2str(eq_style))
                ierr = 1
                CHCKERR(err_msg)
        end select
    end function calc_ang_grid_eq
    
    ! Calculate grid that follows magnetic field lines.
    integer function calc_ang_grid_eq_B(grid_eq,eq) result(ierr)
        use num_vars, only: use_pol_flux_F, use_pol_flux_E, &
            &eq_style, tol_NR
        use grid_vars, only: min_par_X, max_par_X, alpha
        use eq_vars, only: max_flux_p_E, max_flux_t_E
        
        character(*), parameter :: rout_name = 'calc_ang_grid_eq_B'
        
        ! input / output
        type(grid_type), intent(inout) :: grid_eq                               ! equilibrium grid of which to calculate angular part
        type(eq_type), intent(in) :: eq                                         ! equilibrium containing the angular grid
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        real(dp), allocatable :: r_E_loc(:)                                     ! flux in Equilibrium coords.
        real(dp), pointer :: flux_F(:) => null()                                ! flux that the F uses as normal coord.
        real(dp), pointer :: flux_E(:) => null()                                ! flux that the E uses as normal coord.
        real(dp) :: r_F_factor, r_E_factor                                      ! mult. factors for r_F and r_E
        integer :: pmone                                                        ! plus or minus one
        integer :: kd                                                           ! counter
        
        ! initialize ierr
        ierr = 0
        
        ! set up flux_E and plus minus one
        ! Note: this routine is similar  to calc_loc_r, but that routine
        ! cannot  be  used here  because  it  is  possible that  the  FD
        ! quantities are not yet defined.
        ! choose which equilibrium style is being used:
        !   1:  VMEC
        !   2:  HELENA
        select case (eq_style)
            case (1)                                                            ! VMEC
                if (use_pol_flux_E) then                                        ! normalized poloidal flux
                    flux_E => eq%flux_p_E(:,0)
                    r_E_factor = max_flux_p_E
                else                                                            ! normalized toroidal flux
                    flux_E => eq%flux_t_E(:,0)
                    r_E_factor = max_flux_t_E
                end if
                pmone = -1                                                      ! conversion VMEC LH -> RH coord. system
            case (2)                                                            ! HELENA
                flux_E => eq%flux_p_E(:,0)
                r_E_factor = 2*pi
                pmone = 1
            case default
                err_msg = 'No equilibrium style associated with '//&
                    &trim(i2str(eq_style))
                ierr = 1
                CHCKERR(err_msg)
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
        if (use_pol_flux_F) then                                                ! parallel angle theta
            ierr = calc_eqd_grid(grid_eq%theta_F,min_par_X*pi,&
                &max_par_X*pi,1,excl_last=.true.)                               ! first index corresponds to parallel angle
            CHCKERR('')
            do kd = 1,grid_eq%loc_n_r
                grid_eq%zeta_F(:,:,kd) = pmone*eq%q_saf_E(kd,0)*&
                    &grid_eq%theta_F(:,:,kd)
            end do
            grid_eq%zeta_F = grid_eq%zeta_F + alpha
        else                                                                    ! parallel angle zeta
            ierr = calc_eqd_grid(grid_eq%zeta_F,min_par_X*pi,&
                &max_par_X*pi,1,excl_last=.true.)                               ! first index corresponds to parallel angle
            CHCKERR('')
            do kd = 1,grid_eq%loc_n_r
                grid_eq%theta_F(:,:,kd) = pmone*eq%rot_t_E(kd,0)*&
                    &grid_eq%zeta_F(:,:,kd)
            end do
            grid_eq%theta_F = grid_eq%theta_F - alpha
        end if
        
        ! allocate local r_E
        allocate(r_E_loc(size(flux_F)))
        
        ! convert  Flux  coordinates  to  Equilibrium  coordinates  (use
        ! custom flux_E and  flux_F because the Flux  quantities are not
        ! yet calculated)
        ierr = coord_F2E(eq,grid_eq,&
            &grid_eq%loc_r_F,grid_eq%theta_F,grid_eq%zeta_F,&
            &r_E_loc,grid_eq%theta_E,grid_eq%zeta_E,&
            &r_F_array=flux_F/r_F_factor,r_E_array=flux_E/r_E_factor)
        CHCKERR('')
        
        ! test whether r_E_loc indeed corresponds to loc_r_E of eq grid
        if (maxval(abs(grid_eq%loc_r_E-r_E_loc)).gt.10*tol_NR) then
            ierr = 1
            err_msg = 'loc_r_E of equilibrium grid is not recovered'
            CHCKERR(err_msg)
        end if
        
        ! deallocate local variables
        deallocate(r_E_loc)
        nullify(flux_F,flux_E)
        
#if ldebug
        if (debug_calc_ang_grid_eq_B) then
            call writo('Plotting theta_E, theta_F, zeta_E and zeta_F')
            call plot_HDF5('TEST_theta_E','TEST_theta_E',grid_eq%theta_E)
            call plot_HDF5('TEST_theta_F','TEST_theta_F',grid_eq%theta_F)
            call plot_HDF5('TEST_zeta_E','TEST_zeta_E',grid_eq%zeta_E)
            call plot_HDF5('TEST_zeta_F','TEST_zeta_F',grid_eq%zeta_F)
        end if
#endif
    end function calc_ang_grid_eq_B
    
    ! Converts  Flux  coordinates  (r,theta,zeta)_F to  Equilibrium  coordinates
    ! (r,theta,zeta)_E. Optionally,  two arrays  r_F_array and r_E_array  can be
    ! provided,  that define  the mapping  between the  both coordinate  system.
    ! Standard, for E, the poloidal or  toroidal normalized flux is used and for
    ! F, the poloidal or toroidal flux in F coordinates, divided by 2pi.
    integer function coord_F2E_rtz(eq,grid_eq,r_F,theta_F,zeta_F,r_E,&
        &theta_E,zeta_E,r_F_array,r_E_array) result(ierr)                       ! version with r, theta and zeta
        use num_vars, only: tol_NR, eq_style
        use fourier, only: fourier2real, calc_trigon_factors
        use VMEC, only: L_V_c, L_V_s
        use utilities, only: interp_fun
        
        character(*), parameter :: rout_name = 'coord_F2E_rtz'
        
        ! input / output
        type(eq_type), intent(in) :: eq                                         ! equilibrium in which to convert variables
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid (for normal local limits)
        real(dp), intent(in) :: r_F(:), theta_F(:,:,:), zeta_F(:,:,:)           ! Flux coords.
        real(dp), intent(inout) :: r_E(:), theta_E(:,:,:), zeta_E(:,:,:)        ! Equilibrium coords.
        real(dp), intent(in), optional, target :: r_F_array(:), r_E_array(:)    ! optional arrays that define mapping between two coord. systems
        
        ! local variables (also used in child functions)
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer :: n_r, n_par, n_geo                                            ! dimensions of the grid
        real(dp), allocatable :: L_V_c_loc(:,:,:)                               ! local version of L_V_c
        real(dp), allocatable :: L_V_s_loc(:,:,:)                               ! local version of L_V_s
        real(dp) :: theta_V_loc(1,1,1), zeta_V_loc(1,1,1)                       ! local version of theta_V, zeta_V (needs to be 3D matrix)
        real(dp) :: lam(1,1,1)                                                  ! lambda (needs to be 3D matrix)
        real(dp) :: dlam(1,1,1)                                                 ! angular derivative of lambda (needs to be 3D matrix)
        integer :: id, jd, kd                                                   ! counters
        
        ! initialize ierr
        ierr = 0
        
        ! set up array sizes
        n_par = size(theta_E,1)
        n_geo = size(theta_E,2)
        n_r = size(theta_E,3)
        
        ! tests
        if (size(theta_F,1).ne.size(zeta_F,1) .or. &
            &size(theta_F,2).ne.size(zeta_F,2) .or. &
            &size(theta_F,3).ne.size(zeta_F,3) .or. &
            &size(theta_F,3).ne.size(r_F)) then
            ierr = 1
            err_msg = 'theta_F, zeta_F and r_F need to have the correct &
                &dimensions'
            CHCKERR(err_msg)
        end if
        if (n_par.ne.size(zeta_E,1) .or. n_geo.ne.size(zeta_E,2) .or. &
            &n_r.ne.size(zeta_E,3) .or. n_r.ne.size(r_E)) then
            ierr = 1
            err_msg = 'theta_E, zeta_E and r_E need to have the correct &
                &dimensions'
            CHCKERR(err_msg)
        end if
        if (present(r_F_array).neqv.present(r_E_array)) then
            ierr = 1
            err_msg = 'both r_F_array and r_E_array have to be provided'
            CHCKERR(err_msg)
        end if
        
        ! convert normal position
        ierr = coord_F2E_r(grid_eq,eq,r_F,r_E,r_F_array,r_E_array)
        CHCKERR('')
        
        ! choose which equilibrium style is being used:
        !   1:  VMEC
        !   2:  HELENA
        select case (eq_style)
            case (1)                                                            ! VMEC
                ierr = coord_F2E_VMEC()
                CHCKERR('')
            case (2)                                                            ! HELENA
                ! trivial HELENA uses flux coordinates
                theta_E = theta_F
                zeta_E = zeta_F
            case default
                err_msg = 'No equilibrium style associated with '//&
                    &trim(i2str(eq_style))
                ierr = 1
                CHCKERR(err_msg)
        end select
    contains
        integer function coord_F2E_VMEC() result(ierr)
            use utilities, only: calc_zero_NR
            use VMEC, only: mpol, ntor
            
            character(*), parameter :: rout_name = 'coord_F2E_VMEC'
            
            ! local variables
            real(dp), allocatable :: L_V_c_loc_ind(:,:), L_V_s_loc_ind(:,:)     ! individual versions of L_V_c_loc and L_V_s_loc
            real(dp), allocatable :: loc_r_F(:)                                 ! loc_r in F coords.
            
            ! initialize ierr
            ierr = 0
            
            ! set up loc_r_F
            if (present(r_F_array)) then
                allocate(loc_r_F(grid_eq%loc_n_r))
                loc_r_F = r_F_array
            else
                ierr = calc_loc_r(grid_eq,eq,loc_r_F,style=2)
                CHCKERR('')
            end if
            
            ! the toroidal angle is trivial
            zeta_E = - zeta_F                                                   ! conversion VMEC LH -> RH coord. system
            
            ! allocate local copies of L_V_c and L_V_s
            allocate(L_V_c_loc(0:mpol-1,-ntor:ntor,1:1))
            allocate(L_V_s_loc(0:mpol-1,-ntor:ntor,1:1))
            allocate(L_V_c_loc_ind(0:mpol-1,-ntor:ntor))
            allocate(L_V_s_loc_ind(0:mpol-1,-ntor:ntor))
            
            ! loop over all normal points
            do kd = 1,n_r
                ! interpolate L_V_c and L_V_s at requested normal point r_F
                ierr = interp_fun(L_V_c_loc_ind,&
                    &L_V_c(:,:,grid_eq%i_min:grid_eq%i_max,0),r_F(kd),loc_r_F)
                CHCKERR('')
                ierr = interp_fun(L_V_s_loc_ind,&
                    &L_V_s(:,:,grid_eq%i_min:grid_eq%i_max,0),r_F(kd),loc_r_F)
                CHCKERR('')
                
                ! copy individual to array version
                L_V_c_loc(:,:,1) = L_V_c_loc_ind
                L_V_s_loc(:,:,1) = L_V_s_loc_ind
                
                ! the poloidal angle has to be found as the zero of
                !   f = theta_F - theta_E - lambda
                ! loop over all angular points
                do jd = 1,n_geo
                    do id = 1,n_par
                        ! calculate zero of f
                        ierr = calc_zero_NR(theta_E(id,jd,kd),fun_pol,&
                            &dfun_pol,theta_F(id,jd,kd))                        ! use theta_F as guess for theta_E
                        CHCKERR('')
                        
                        ! do a check whether the result is indeed zero
                        if (abs(fun_pol(theta_E(id,jd,kd))).gt.tol_NR*100) then
                            err_msg = 'In coord_F2E_VMEC, calculating f as a &
                                &check, using the theta_E that is the solution&
                                & of f = 0, yields a calculated f that &
                                &deviates from 0 by '//trim(r2strt(&
                                &100*abs(fun_pol(theta_E(id,jd,kd)))))//'%'
                            ierr = 1
                            CHCKERR(err_msg)
                        end if
                    end do
                end do
            end do
            
            ! deallocate variables
            deallocate(L_V_c_loc,L_V_s_loc)
            deallocate(L_V_c_loc_ind,L_V_s_loc_ind)
        end function coord_F2E_VMEC
        
        ! function that returns f = theta_F  - theta_V - lambda. It uses theta_F
        ! and  zeta_E (=  zeta_V)  from the  parent function,  lam  and dlam  to
        ! contain the variable  lambda or its derivative, as  well as L_V_s_loc,
        ! L_V_c_loc, theta_V_loc, zeta_V_loc, id, jd and kd.
        function fun_pol(theta_E_in)
            character(*), parameter :: rout_name = 'fun_pol'
            
            ! input / output
            real(dp) :: fun_pol
            real(dp), intent(in) :: theta_E_in
            
            ! local variables
            real(dp), allocatable :: trigon_factors_loc(:,:,:,:,:,:)            ! trigonometric factor cosine for the inverse fourier transf.
            
            ! initialize fun_pol
            fun_pol = 0.0_dp
            
            ! set up local copies of theta_V_loc, zeta_V_loc
            theta_V_loc = theta_E_in
            zeta_V_loc = zeta_E(id,jd,kd)
            
            ! transform lambda from Fourier space to real space
            ! calculate the (co)sines
            ierr = calc_trigon_factors(theta_V_loc,zeta_V_loc,&
                &trigon_factors_loc)
            CHCKERR('')
            
            ! calculate lambda
            ierr = fourier2real(L_V_c_loc,L_V_s_loc,trigon_factors_loc,lam)
            CHCKERR('')
            
            ! calculate the output function
            fun_pol = theta_F(id,jd,kd) - theta_E_in - lam(1,1,1)
            
            ! deallocate trigoniometric factors
            deallocate(trigon_factors_loc)
        end function fun_pol
        
        ! function that  returns df/dtheta_V  = -1  - dlambda/dtheta_V.  It uses
        ! theta_F and zeta_E  (= zeta_V) from the parent function,  lam and dlam
        ! to  contain  the  variable  lambda  or  its  derivative,  as  well  as
        ! L_V_s_loc, L_V_c_loc, theta_V_loc, zeta_V_loc, id, jd and kd.
        function dfun_pol(theta_E_in)
            character(*), parameter :: rout_name = 'dfun_pol'
            
            ! input / output
            real(dp) :: dfun_pol
            real(dp), intent(in) :: theta_E_in
            
            ! local variables
            real(dp), allocatable :: trigon_factors_loc(:,:,:,:,:,:)            ! trigonometric factor cosine for the inverse fourier transf.
            
            ! initialize dfun_pol
            dfun_pol = 0.0_dp
            
            ! set up local copies of theta_V_loc, zeta_V_loc
            theta_V_loc = theta_E_in
            zeta_V_loc = zeta_E(id,jd,kd)
            
            ! transform lambda from Fourier space to real space
            ! calculate the (co)sines
            ierr = calc_trigon_factors(theta_V_loc,zeta_V_loc,&
                &trigon_factors_loc)
            CHCKERR('')
            
            ! calculate lambda
            ierr = fourier2real(L_V_c_loc,L_V_s_loc,trigon_factors_loc,dlam,[1,0])
            CHCKERR('')
            
            ! calculate the output function
            dfun_pol = -1._dp - dlam(1,1,1)
            
            ! deallocate trigoniometric factors
            deallocate(trigon_factors_loc)
        end function dfun_pol
    end function coord_F2E_rtz
    integer function coord_F2E_r(grid_eq,eq,r_F,r_E,r_F_array,r_E_array) &
        &result(ierr)                                                           ! version with only r
        use utilities, only: interp_fun
        
        character(*), parameter :: rout_name = 'coord_F2E_r'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid (for normal local limits)
        type(eq_type), intent(in) :: eq                                         ! equilibrium in which to convert variables
        real(dp), intent(in) :: r_F(:)                                          ! Flux coords.
        real(dp), intent(inout) :: r_E(:)                                       ! Equilibrium coords.
        real(dp), intent(in), optional, target :: r_F_array(:), r_E_array(:)    ! optional arrays that define mapping between two coord. systems
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer :: n_r                                                          ! dimension of the grid
        integer :: kd                                                           ! counter
        real(dp), allocatable :: loc_r_E(:), loc_r_F(:)                         ! loc_r in E and F coords.
        
        ! initialize ierr
        ierr = 0
        
        ! set up array sizes
        n_r = size(r_F)
        
        ! tests
        if (n_r.ne.size(r_E)) then
            ierr = 1
            err_msg = 'r_F and r_E need to have the correct dimensions'
            CHCKERR(err_msg)
        end if
        if (present(r_F_array).neqv.present(r_E_array)) then
            ierr = 1
            err_msg = 'both r_F_array and r_E_array have to be provided'
            CHCKERR(err_msg)
        end if
        
        ! set up loc_r_E and loc_r_F
        if (present(r_F_array)) then
            allocate(loc_r_E(grid_eq%loc_n_r))
            allocate(loc_r_F(grid_eq%loc_n_r))
            loc_r_E = r_E_array
            loc_r_F = r_F_array
        else
            ierr = calc_loc_r(grid_eq,eq,loc_r_E,style=1)
            CHCKERR('')
            ierr = calc_loc_r(grid_eq,eq,loc_r_F,style=2)
            CHCKERR('')
        end if
        
        ! convert normal position
        do kd = 1,n_r
            ierr = interp_fun(r_E(kd),loc_r_E,r_F(kd),loc_r_F)
            CHCKERR('')
            r_E(kd) = r_E(kd)
        end do
    end function coord_F2E_r
    
    ! Converts  Equilibrium  coordinates  (r,theta,zeta)_E to  Flux  coordinates
    ! (r,theta,zeta)_F. Optionally,  two arrays  r_E_array and r_F_array  can be
    ! provided,  that define  the mapping  between the  both coordinate  system.
    ! Standard, for E, the poloidal or  toroidal normalized flux is used and for
    ! F, the poloidal or toroidal flux in F coordinates, divided by 2pi.
    integer function coord_E2F_rtz(eq,grid_eq,r_E,theta_E,zeta_E,r_F,&
        &theta_F,zeta_F,r_E_array,r_F_array) result(ierr)                       ! version with r, theta and zeta
        use num_vars, only: eq_style
        use utilities, only: interp_fun
        
        character(*), parameter :: rout_name = 'coord_E2F_rtz'
        
        ! input / output
        type(eq_type), intent(in) :: eq                                         ! equilibrium in which to convert variables
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid (for normal local limits)
        real(dp), intent(in) :: r_E(:), theta_E(:,:,:), zeta_E(:,:,:)           ! Equilibrium coords.
        real(dp), intent(inout) :: r_F(:), theta_F(:,:,:), zeta_F(:,:,:)        ! Flux coords.
        real(dp), intent(in), optional, target :: r_E_array(:), r_F_array(:)    ! optional arrays that define mapping between two coord. systems
        
        ! local variables (also used in child functions)
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer :: n_r, n_par, n_geo                                            ! dimensions of the grid
        integer :: kd                                                           ! counter
        
        ! initialize ierr
        ierr = 0
        
        ! set up array sizes
        n_par = size(theta_E,1)
        n_geo = size(theta_E,2)
        n_r = size(theta_E,3)
        
        ! tests
        if (size(theta_F,1).ne.size(zeta_F,1) .or. &
            &size(theta_F,2).ne.size(zeta_F,2) .or. &
            &size(theta_F,3).ne.size(zeta_F,3) .or. &
            &size(theta_F,3).ne.size(r_F)) then
            ierr = 1
            err_msg = 'theta_F, zeta_F and r_F need to have the correct &
                &dimensions'
            CHCKERR(err_msg)
        end if
        if (n_par.ne.size(zeta_E,1) .or. n_geo.ne.size(zeta_E,2) .or. &
            &n_r.ne.size(zeta_E,3) .or. n_r.ne.size(r_E)) then
            ierr = 1
            err_msg = 'theta_E, zeta_E and r_E need to have the correct &
                &dimensions'
            CHCKERR(err_msg)
        end if
        if (present(r_F_array).neqv.present(r_E_array)) then
            ierr = 1
            err_msg = 'both r_F_array and r_E_array have to be provided'
            CHCKERR(err_msg)
        end if
        
        ! convert normal position
        ierr = coord_E2F_r(grid_eq,eq,r_E,r_F,r_E_array,r_F_array)
        CHCKERR('')
        
        ! choose which equilibrium style is being used:
        !   1:  VMEC
        !   2:  HELENA
        select case (eq_style)
            case (1)                                                            ! VMEC
                ierr = coord_E2F_VMEC()
                CHCKERR('')
            case (2)                                                            ! HELENA
                ! trivial HELENA uses flux coordinates
                theta_F = theta_E
                zeta_F = zeta_E
            case default
                err_msg = 'No equilibrium style associated with '//&
                    &trim(i2str(eq_style))
                ierr = 1
                CHCKERR(err_msg)
        end select
    contains
        integer function coord_E2F_VMEC() result(ierr)
            use VMEC, only: mpol, ntor, L_V_c, L_V_s
            use fourier, only: calc_trigon_factors, fourier2real
            
            character(*), parameter :: rout_name = 'coord_E2F_VMEC'
            
            ! local variables
            real(dp), allocatable :: L_V_c_loc(:,:,:), L_V_s_loc(:,:,:)         ! local version of L_V_c and L_V_s
            real(dp), allocatable :: trigon_factors_loc(:,:,:,:,:,:)            ! trigonometric factor cosine for the inverse fourier transf.
            real(dp), allocatable :: lam(:,:,:)                                 ! lambda
            real(dp), allocatable :: loc_r_E(:)                                 ! loc_r in E coords.
            
            ! initialize ierr
            ierr = 0
            
            ! set up loc_r_E
            if (present(r_E_array)) then
                allocate(loc_r_E(grid_eq%loc_n_r))
                loc_r_E = r_E_array
            else
                ierr = calc_loc_r(grid_eq,eq,loc_r_E,style=1)
                CHCKERR('')
            end if
            
            ! the toroidal angle is trivial
            zeta_F = - zeta_E                                                   ! conversion VMEC LH -> RH coord. system
            
            ! allocate local copies of L_V_c and L_V_s and lambda
            allocate(L_V_c_loc(0:mpol-1,-ntor:ntor,1:n_r))
            allocate(L_V_s_loc(0:mpol-1,-ntor:ntor,1:n_r))
            allocate(lam(n_par,n_geo,n_r))
            
            ! interpolate L_V_c and L_V_s at requested normal point r_E
            ! loop over all normal points
            do kd = 1,n_r
                ierr = interp_fun(L_V_c_loc(:,:,kd),&
                    &L_V_c(:,:,grid_eq%i_min:grid_eq%i_max,0),r_E(kd),loc_r_E)
                CHCKERR('')
                ierr = interp_fun(L_V_s_loc(:,:,kd),&
                    &L_V_s(:,:,grid_eq%i_min:grid_eq%i_max,0),r_E(kd),loc_r_E)
                CHCKERR('')
            end do
            
            ! calculate the (co)sines to transform lambda to real space
            ierr = calc_trigon_factors(theta_E,zeta_E,trigon_factors_loc)
            CHCKERR('')
            
            ! calculate lambda
            ierr = fourier2real(L_V_c_loc,L_V_s_loc,trigon_factors_loc,lam)
            CHCKERR('')
            
            ! the poloidal angle has to be found as
            !   theta_F = theta_E + lambda
            theta_F = theta_E + lam
            
            ! deallocate variables
            deallocate(trigon_factors_loc)
            deallocate(L_V_c_loc,L_V_s_loc)
            deallocate(lam)
        end function coord_E2F_VMEC
    end function coord_E2F_rtz
    integer function coord_E2F_r(grid_eq,eq,r_E,r_F,r_E_array,r_F_array) &
        &result(ierr)                                                           ! version with only r
        use utilities, only: interp_fun
        
        character(*), parameter :: rout_name = 'coord_E2F_r'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid (for normal local limits)
        type(eq_type), intent(in) :: eq                                         ! equilibrium in which to convert variables
        real(dp), intent(in) :: r_E(:)                                          ! Equilibrium coords.
        real(dp), intent(inout) :: r_F(:)                                       ! Flux coords.
        real(dp), intent(in), optional, target :: r_E_array(:), r_F_array(:)    ! optional arrays that define mapping between two coord. systems
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer :: n_r                                                          ! dimension of the grid
        integer :: kd                                                           ! counter
        real(dp), allocatable :: loc_r_E(:), loc_r_F(:)                         ! loc_r in E and F coords.
        
        ! initialize ierr
        ierr = 0
        
        ! set up array sizes
        n_r = size(r_E)
        
        ! tests
        if (n_r.ne.size(r_F)) then
            ierr = 1
            err_msg = 'r_E and r_F need to have the correct dimensions'
            CHCKERR(err_msg)
        end if
        if (present(r_E_array).neqv.present(r_F_array)) then
            ierr = 1
            err_msg = 'both r_E_array and r_F_array have to be provided'
            CHCKERR(err_msg)
        end if
        
        ! set up loc_r_E and loc_r_F
        if (present(r_F_array)) then
            allocate(loc_r_E(grid_eq%loc_n_r))
            allocate(loc_r_F(grid_eq%loc_n_r))
            loc_r_E = r_E_array
            loc_r_F = r_F_array
        else
            ierr = calc_loc_r(grid_eq,eq,loc_r_E,style=1)
            CHCKERR('')
            ierr = calc_loc_r(grid_eq,eq,loc_r_F,style=2)
            CHCKERR('')
        end if
        
        ! convert normal position
        do kd = 1,n_r
            ierr = interp_fun(r_F(kd),loc_r_F,r_E(kd),loc_r_E)
            CHCKERR('')
        end do
    end function coord_E2F_r
    
    ! Calculates X,Y  and Z on a  grid, which should have  the local equilibrium
    ! grid angles  set up in E(quilibrium)  coordinates. The total r_E,  and the
    ! F(lux) variables are ignored.
    ! If VMEC is the equilibrium  model, this routine also optionally calculates
    ! lambda on the grid, as this is  also needed some times. If HELENA is used,
    ! this variable is not used.
    ! X, Y, Z and  optionally L need to have the correct  dimensions if they are
    ! allocated.  They  can  also  be  passed unallocated,  in  which  case  the
    ! allocation happens automatically.
    integer function calc_XYZ_grid(grid,X,Y,Z,L) result(ierr)
        use num_vars, only: eq_style, use_normalization
        use utilities, only: interp_fun, round_with_tol
        use eq_vars, only: R_0
        
        character(*), parameter :: rout_name = 'calc_XYZ_grid'
        
        ! input / output
        type(grid_type), intent(in) :: grid                                     ! grid for which to calculate X, Y, Z and optionally L
        real(dp), intent(inout) :: X(:,:,:), Y(:,:,:), Z(:,:,:)                 ! X, Y and Z of grid
        real(dp), intent(inout), optional :: L(:,:,:)                           ! lambda of grid
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! test
        if (size(X,1).ne.grid%n(1) .or. size(X,2).ne.grid%n(2) .or. &
            &size(X,3).ne.grid%loc_n_r) then
            ierr = 1
            err_msg =  'X needs to have the correct dimensions'
            CHCKERR(err_msg)
        end if
        if (size(Y,1).ne.grid%n(1) .or. size(Y,2).ne.grid%n(2) .or. &
            &size(Y,3).ne.grid%loc_n_r) then
            ierr = 1
            err_msg =  'Y needs to have the correct dimensions'
            CHCKERR(err_msg)
        end if
        if (size(Z,1).ne.grid%n(1) .or. size(Z,2).ne.grid%n(2) .or. &
            &size(Z,3).ne.grid%loc_n_r) then
            ierr = 1
            err_msg =  'Z needs to have the correct dimensions'
            CHCKERR(err_msg)
        end if
        if (present(L)) then
            if (size(L,1).ne.grid%n(1) .or. size(L,2).ne.grid%n(2) .or. &
                &size(L,3).ne.grid%loc_n_r) then
                ierr = 1
                err_msg =  'L needs to have the correct dimensions'
                CHCKERR(err_msg)
            end if
        end if
        
        ! choose which equilibrium style is being used:
        !   1:  VMEC
        !   2:  HELENA
        select case (eq_style)
            case (1)                                                            ! VMEC
                ierr = calc_XYZ_grid_VMEC(grid,X,Y,Z,L)
                CHCKERR('')
            case (2)                                                            ! HELENA
                ierr = calc_XYZ_grid_HEL(grid,X,Y,Z)
                CHCKERR('')
            case default
                err_msg = 'No equilibrium style associated with '//&
                    &trim(i2str(eq_style))
                ierr = 1
                CHCKERR(err_msg)
        end select
        
        ! if normalized, return to physical value
        if (use_normalization) then
            X = X*R_0
            Y = Y*R_0
            Z = Z*R_0
        end if
    contains
        ! VMEC version
        integer function calc_XYZ_grid_VMEC(grid,X,Y,Z,L) result(ierr)
            use VMEC, only: R_V_c, R_V_s, Z_V_c, Z_V_s, L_V_c, L_V_s, mpol, ntor
            use fourier, only: calc_trigon_factors, fourier2real
            
            character(*), parameter :: rout_name = 'calc_XYZ_grid_VMEC'
            
            ! input / output
            type(grid_type), intent(in) :: grid                                 ! grid for which to calculate X, Y, Z and optionally L
            real(dp), intent(inout) :: X(:,:,:), Y(:,:,:), Z(:,:,:)             ! X, Y and Z of grid
            real(dp), intent(inout), optional :: L(:,:,:)                       ! lambda of grid
            
            ! local variables
            integer :: kd, m, n                                                 ! counters
            real(dp), allocatable :: R_V_c_int(:,:,:), R_V_s_int(:,:,:)         ! interpolated version of R_V_c and R_V_s
            real(dp), allocatable :: Z_V_c_int(:,:,:), Z_V_s_int(:,:,:)         ! interpolated version of Z_V_c and Z_V_s
            real(dp), allocatable :: L_V_c_int(:,:,:), L_V_s_int(:,:,:)         ! interpolated version of L_V_c and L_V_s
            real(dp), allocatable :: trigon_factors(:,:,:,:,:,:)                ! trigonometric factor cosine for the inverse fourier transf.
            real(dp), allocatable :: R(:,:,:)                                   ! R in Cylindrical coordinates
            
            ! initialize ierr
            ierr = 0
            
            ! set up interpolated R_V_c_int, ..
            allocate(R_V_c_int(0:mpol-1,-ntor:ntor,1:grid%loc_n_r))
            allocate(R_V_s_int(0:mpol-1,-ntor:ntor,1:grid%loc_n_r))
            allocate(Z_V_c_int(0:mpol-1,-ntor:ntor,1:grid%loc_n_r))
            allocate(Z_V_s_int(0:mpol-1,-ntor:ntor,1:grid%loc_n_r))
            if (present(L)) then
                allocate(L_V_c_int(0:mpol-1,-ntor:ntor,1:grid%loc_n_r))
                allocate(L_V_s_int(0:mpol-1,-ntor:ntor,1:grid%loc_n_r))
            end if
            
            ! interpolate VMEC tables for every requested normal point
            do kd = 1,grid%loc_n_r                                              ! loop over all local normal points
                ! interpolate the VMEC tables in normal direction
                ! Note: VMEC uses an equidistant grid normalized from 0 to 1
                do n = -ntor,ntor
                    do m = 0,mpol-1
                        ierr = interp_fun(R_V_c_int(m,n,kd),R_V_c(m,n,:,0),&
                            &grid%loc_r_E(kd))
                        CHCKERR('')
                        ierr = interp_fun(R_V_s_int(m,n,kd),R_V_s(m,n,:,0),&
                            &grid%loc_r_E(kd))
                        CHCKERR('')
                        ierr = interp_fun(Z_V_c_int(m,n,kd),Z_V_c(m,n,:,0),&
                            &grid%loc_r_E(kd))
                        CHCKERR('')
                        ierr = interp_fun(Z_V_s_int(m,n,kd),Z_V_s(m,n,:,0),&
                            &grid%loc_r_E(kd))
                        CHCKERR('')
                        if (present(L)) then
                            ierr = interp_fun(L_V_c_int(m,n,kd),L_V_c(m,n,:,0),&
                                &grid%loc_r_E(kd))
                            CHCKERR('')
                            ierr = interp_fun(L_V_s_int(m,n,kd),L_V_s(m,n,:,0),&
                                &grid%loc_r_E(kd))
                            CHCKERR('')
                        end if
                    end do
                end do
            end do
            
            ! calculate trigonometric factors
            ierr = calc_trigon_factors(grid%theta_E,grid%zeta_E,trigon_factors)
            CHCKERR('')
            
            ! allocate R
            allocate(R(grid%n(1),grid%n(2),grid%loc_n_r))
            
            ! inverse fourier transform with trigonometric factors
            ierr = fourier2real(R_V_c_int,R_V_s_int,trigon_factors,R,[0,0])
            CHCKERR('')
            ierr = fourier2real(Z_V_c_int,Z_V_s_int,trigon_factors,Z,[0,0])
            CHCKERR('')
            if (present(L)) then
                ierr = fourier2real(L_V_c_int,L_V_s_int,trigon_factors,L,[0,0])
                CHCKERR('')
            end if
            
            ! deallocate
            deallocate(trigon_factors)
            
            ! transform cylindrical to cartesian
            ! (the geometrical zeta is the inverse of VMEC zeta)
            X = R*cos(-grid%zeta_E)
            Y = R*sin(-grid%zeta_E)
            
            ! deallocate
            deallocate(R_V_c_int,R_V_s_int,Z_V_c_int,Z_V_s_int)
            if (present(L)) deallocate(L_V_c_int,L_V_s_int)
            deallocate(R)
        end function calc_XYZ_grid_VMEC
        
        ! HELENA version
        integer function calc_XYZ_grid_HEL(grid,X,Y,Z) result(ierr)
            use HELENA, only: R_H, Z_H, chi_H, ias, flux_p_H
            
            character(*), parameter :: rout_name = 'calc_XYZ_grid_HEL'
            
            ! input / output
            type(grid_type), intent(in) :: grid                                 ! grid for which to calculate X, Y, Z and optionally L
            real(dp), intent(inout) :: X(:,:,:), Y(:,:,:), Z(:,:,:)             ! X, Y and Z of grid
            
            ! local variables
            integer :: id, jd, kd                                               ! counters
            real(dp), allocatable :: R_H_int(:), Z_H_int(:)                     ! R and Z at interpolated normal value
            real(dp), allocatable :: R(:,:,:)                                   ! R in Cylindrical coordinates
            real(dp) :: theta_loc                                               ! local copy of theta_E
            
            ! initialize ierr
            ierr = 0
            
            ! allocate R
            allocate(R(grid%n(1),grid%n(2),grid%loc_n_r))
            
            ! set up interpolated R and Z
            allocate(R_H_int(size(R_H,1)),Z_H_int(size(Z_H,1)))
            
            ! interpolate HELENA output  R_H and Z_H for  every requested normal
            ! point
            ! Note:  R_H and  Z_H  are not  adapted to  the  parallel grid,  but
            ! tabulated in the original HELENA poloidal grid.
            ! Note:  HELENA uses poloidal  flux as normal coordinate  divided by
            ! 2pi
            ! interpolate the HELENA tables in normal direction
            do kd = 1,grid%loc_n_r                                              ! loop over all normal points
                do id = 1,size(R_H,1)
                    ierr = interp_fun(R_H_int(id),R_H(id,:),grid%loc_r_E(kd),&
                        &x=flux_p_H/(2*pi))
                    CHCKERR('')
                end do
                do id = 1,size(Z_H,1)
                    ierr = interp_fun(Z_H_int(id),Z_H(id,:),grid%loc_r_E(kd),&
                        &x=flux_p_H/(2*pi))
                    CHCKERR('')
                end do
                ! loop over toroidal points
                do jd = 1,grid%n(2)
                    ! interpolate at the requested poloidal points
                    do id = 1,grid%n(1)
                        theta_loc = grid%theta_E(id,jd,kd)
                        ! add or subtract 2pi to  the parallel angle until it is
                        ! at least 0 to get principal range 0..2pi
                        if (theta_loc.lt.0._dp) then
                            do while (theta_loc.lt.0._dp)
                                theta_loc = theta_loc + 2*pi
                            end do
                        else if (theta_loc.gt.2*pi) then
                            do while (theta_loc.gt.2._dp*pi)
                                theta_loc = theta_loc - 2*pi
                            end do
                        end if
                        ! Interpolate  the HELENA  variables poloidally,  taking
                        ! into account the possible symmetry
                        if (ias.eq.0 .and. theta_loc.gt.pi) then
                            ierr = interp_fun(R(id,jd,kd),R_H_int,&
                                &2*pi-theta_loc,x=chi_H)
                            CHCKERR('')
                            ierr = interp_fun(Z(id,jd,kd),-Z_H_int,&
                                &2*pi-theta_loc,x=chi_H)                        ! Z at second half of period is inverted
                            CHCKERR('')
                        else
                            ierr = interp_fun(R(id,jd,kd),R_H_int,theta_loc,&
                                &x=chi_H)
                            CHCKERR('')
                            ierr = interp_fun(Z(id,jd,kd),Z_H_int,theta_loc,&
                                &x=chi_H)
                            CHCKERR('')
                        end if
                    end do
                end do
            end do
            
            ! calculate X and Y, transforming cylindrical to cartesian
            X = R*cos(grid%zeta_E)
            Y = R*sin(grid%zeta_E)
            
            ! deallocate
            deallocate(R)
        end function calc_XYZ_grid_HEL
    end function calc_XYZ_grid

    ! Calculate grid of equidistant points,  where optionally the last point can
    ! be excluded.
    integer function calc_eqd_grid_3D(var,min_grid,max_grid,grid_dim,&
        &excl_last) result(ierr)                                                ! 3D version
        character(*), parameter :: rout_name = 'calc_eqd_grid_3D'
        
        ! input and output
        real(dp), intent(inout) :: var(:,:,:)                                   ! output
        real(dp), intent(in) :: min_grid, max_grid                              ! min. and max. of angles [pi]
        integer, intent(in) :: grid_dim                                         ! in which dimension to create the grid
        logical, intent(in), optional :: excl_last                              ! .true. if last point excluded
        
        ! local variables
        integer :: id                                                           ! counter
        character(len=max_str_ln) :: err_msg                                    ! error message
        logical :: excl_last_loc                                                ! local copy of excl_last
        integer :: grid_size                                                    ! nr. of points
        integer :: grid_size_mod                                                ! grid_size or grid_size + 1
        
        ! initialize ierr
        ierr = 0
        
        ! tests
        if (grid_dim.lt.1 .or. grid_dim.gt.3) then
            ierr = 1
            err_msg = 'grid_dim has to point to a dimension going from 1 to 3'
            CHCKERR(err_msg)
        end if
        
        ! set up grid_size
        grid_size = size(var,grid_dim)
        
        ! test some values
        if (grid_size.lt.1) then
            err_msg = 'The angular array has to have a length of at &
                &least 1'
            ierr = 1
            CHCKERR(err_msg)
        end if
        if (min_grid.gt.max_grid) then
            err_msg = 'The minimum angle has to be smaller than the &
                &maximum angle'
            ierr = 1
            CHCKERR(err_msg)
        end if
        
        ! set up local excl_last
        excl_last_loc = .false.
        if (present(excl_last)) excl_last_loc = excl_last
        
        ! add 1 to modified grid size if last point is to be excluded
        if (excl_last_loc) then
            grid_size_mod = grid_size + 1
        else
            grid_size_mod = grid_size
        end if
        
        ! initialize output vector
        var = 0.0_dp
        
        ! calculate grid points
        if (grid_size.eq.1) then
            var = min_grid                                                      ! = max_grid
        else
            if (grid_dim.eq.1) then
                do id = 1,grid_size
                    var(id,:,:) = min_grid + &
                        &(max_grid-min_grid)*(id-1)/(grid_size_mod-1)
                end do
            else if (grid_dim.eq.2) then
                do id = 1,grid_size
                    var(:,id,:) = min_grid + &
                        &(max_grid-min_grid)*(id-1)/(grid_size_mod-1)
                end do
            else
                do id = 1,grid_size
                    var(:,:,id) = min_grid + &
                        &(max_grid-min_grid)*(id-1)/(grid_size_mod-1)
                end do
            end if
        end if
    end function calc_eqd_grid_3D
    integer function calc_eqd_grid_1D(var,min_grid,max_grid,excl_last) &
        &result(ierr)                                                           ! 1D version
        character(*), parameter :: rout_name = 'calc_eqd_grid_1D'
        
        ! input and output
        real(dp), intent(inout) :: var(:)                                       ! output
        real(dp), intent(in) :: min_grid, max_grid                              ! min. and max. of angles [pi]
        logical, intent(in), optional :: excl_last                              ! .true. if last point excluded
        
        ! local variables
        real(dp), allocatable :: var_3D(:,:,:)                                  ! 3D version of var
        
        ! initialize ierr
        ierr = 0
        
        ! set up var_3D
        allocate(var_3D(size(var),1,1))
        
        ! call 3D version
        ierr = calc_eqd_grid_3D(var_3D,min_grid,max_grid,1,excl_last)
        CHCKERR('')
        
        ! update var
        var = var_3D(:,1,1)
    end function calc_eqd_grid_1D
    
    ! plots the grid in real space
    ! The  equilibrium grid  should contain  the fieldline-oriented  angles with
    ! ang_1 the parallel angle and ang_2 the field line label.
    ! Note:  This  routine  does  not  use  n_theta_plot  and  n_zeta_plot  from
    ! num_vars, but instead  temporarily overwrites them with its  own, since it
    ! has to be 3D also in the axisymmetric case.
    ! [MPI] Collective call
    integer function plot_grid_real(grid) result(ierr)
        use num_vars, only: rank, n_procs, no_plots, n_theta_plot, n_zeta_plot
        use grid_vars, only: create_grid, dealloc_grid
        
        character(*), parameter :: rout_name = 'plot_grid_real'
        
        ! input / output
        type(grid_type), intent(in) :: grid                                     ! fieldline-oriented grid
        
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
        integer :: id, jd                                                       ! counters
        integer :: n_theta_plot_old                                             ! backup of n_theta_plot
        integer :: n_zeta_plot_old                                              ! backup of n_zeta_plot
        character(len=*), parameter :: anim_name = &
            &'Magnetic field in flux surfaces'                                  ! name of animation
        
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
        n_theta_plot = 40
        n_zeta_plot = 160
        
        ! extend grid
        ierr = extend_grid_E(grid,grid_ext)
        CHCKERR('')
        
        ! restore n_theta_plot and n_zeta_plot
        n_theta_plot = n_theta_plot_old
        n_zeta_plot = n_zeta_plot_old
        
        ! trim extended grid into plot grid
        ierr = trim_grid(grid_ext,grid_plot)
        CHCKERR('')
        
        ! 1. plot flux surfaces
        call writo('writing flux surfaces')
        
        ! calculate X_1,Y_1 and Z_1
        allocate(X_1(grid_plot%n(1),grid_plot%n(2),grid_plot%loc_n_r))
        allocate(Y_1(grid_plot%n(1),grid_plot%n(2),grid_plot%loc_n_r))
        allocate(Z_1(grid_plot%n(1),grid_plot%n(2),grid_plot%loc_n_r))
        ierr = calc_XYZ_grid(grid_plot,X_1,Y_1,Z_1)
        CHCKERR('')
        
        ! dealloc grid
        call dealloc_grid(grid_plot)
        
        ! 2. plot field lines
        call writo('writing field lines')
        
        ! trim grid into plot grid
        ierr = trim_grid(grid,grid_plot)
        CHCKERR('')
        
        ! calculate X_2,Y_2 and Z_2
        allocate(X_2(grid_plot%n(1),grid_plot%n(2),grid_plot%loc_n_r))
        allocate(Y_2(grid_plot%n(1),grid_plot%n(2),grid_plot%loc_n_r))
        allocate(Z_2(grid_plot%n(1),grid_plot%n(2),grid_plot%loc_n_r))
        ierr = calc_XYZ_grid(grid_plot,X_2,Y_2,Z_2)
        CHCKERR('')
        
        ! dealloc grids
        call dealloc_grid(grid_plot)
        call dealloc_grid(grid_ext)
        
        ! get pointers to full X, Y and Z
        ! The reason for  this is that the  plot is not as simple  as usual, and
        ! also  efficiency  is not  the  biggest  priority. Therefore,  all  the
        ! plotting of the local is handled by a single process, the master.
        call get_full_XYZ(X_1,Y_1,Z_1,X_1_tot,Y_1_tot,Z_1_tot,'flux surfaces')
        call get_full_XYZ(X_2,Y_2,Z_2,X_2_tot,Y_2_tot,Z_2_tot,'field lines')
        
        ierr = plot_grid_real_HDF5(X_1_tot,X_2_tot,Y_1_tot,Y_2_tot,&
            &Z_1_tot,Z_2_tot,anim_name)
        CHCKERR('')
        
        ! clean up
        nullify(X_1_tot,Y_1_tot,Z_1_tot)
        nullify(X_2_tot,Y_2_tot,Z_2_tot)
        
        call lvl_ud(-1)
        
        call writo('Done plotting magnetic field and flux surfaces')
    contains
        ! get pointer to full plotting variables X, Y and Z
        subroutine get_full_XYZ(X,Y,Z,X_tot,Y_tot,Z_tot,merge_name)
            use MPI_utilities, only: get_ser_var
            
            ! input / output
            real(dp), intent(in), target :: X(:,:,:), Y(:,:,:), Z(:,:,:)        ! X, Y and Z of either flux surfaces or magnetic field lines
            real(dp), intent(inout), pointer :: X_tot(:,:,:)                    ! pointer to full X
            real(dp), intent(inout), pointer :: Y_tot(:,:,:)                    ! pointer to full Y
            real(dp), intent(inout), pointer :: Z_tot(:,:,:)                    ! pointer to full Z
            character(len=*) :: merge_name                                      ! name of variable to be merged
            
            ! local variables
            real(dp), allocatable :: loc_XYZ(:)                                 ! local copy of X, Y or Z in an angular point
            real(dp), allocatable :: ser_loc_XYZ(:)                             ! serial copy of loc_XYZ
            integer, allocatable :: tot_dim(:)                                  ! total dimensions for plot 
            
            ! merge plots for flux surfaces if more than one process
            if (n_procs.gt.1) then                                              ! merge local plots
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
                
                allocate(loc_XYZ(sum(tot_dim)))
                do jd = 1,size(X,2)
                    do id = 1,size(X,1)
                        loc_XYZ = X(id,jd,:)
                        ierr = get_ser_var(loc_XYZ,ser_loc_XYZ)
                        CHCKERR('')
                        if (rank.eq.0) X_tot(id,jd,:) = ser_loc_XYZ
                        loc_XYZ = Y(id,jd,:)
                        ierr = get_ser_var(loc_XYZ,ser_loc_XYZ)
                        CHCKERR('')
                        if (rank.eq.0) Y_tot(id,jd,:) = ser_loc_XYZ
                        loc_XYZ = Z(id,jd,:)
                        ierr = get_ser_var(loc_XYZ,ser_loc_XYZ)
                        CHCKERR('')
                        if (rank.eq.0) Z_tot(id,jd,:) = ser_loc_XYZ
                    end do
                end do
            else                                                                ! just point
                X_tot => X
                Y_tot => Y
                Z_tot => Z
            end if
        end subroutine get_full_XYZ
        
        ! Plot with HDF5
        integer function plot_grid_real_HDF5(X_1,X_2,Y_1,Y_2,Z_1,Z_2,&
            &anim_name) result(ierr)
            use HDF5_ops, only: open_HDF5_file, add_HDF5_item, &
                &print_HDF5_top, print_HDF5_geom, print_HDF5_3D_data_item, &
                &print_HDF5_grid, close_HDF5_file
            use HDF5_vars, only: XML_str_type, HDF5_file_type
            
            character(*), parameter :: rout_name = 'plot_grid_real_HDF5'
            
            ! input / output
            real(dp), intent(in) :: X_1(:,:,:), Y_1(:,:,:), Z_1(:,:,:)          ! X, Y and Z of surface in Axisymmetric coordinates
            real(dp), intent(in) :: X_2(:,:,:), Y_2(:,:,:), Z_2(:,:,:)          ! X, Y and Z of magnetic field lines in Axisymmetric coordinates
            character(len=*), intent(in) :: anim_name                           ! name of animation
            
            ! local variables
            character(len=max_str_ln) :: plot_title(2)                          ! plot title for flux surface and for field lines
            integer :: id                                                       ! counter
            type(HDF5_file_type) :: file_info                                   ! file info
            type(XML_str_type) :: grids(2)                                      ! grid in spatial collection (flux surface, field line)
            type(XML_str_type), allocatable :: space_col_grids(:)               ! grid with space collection at different times
            type(XML_str_type) :: time_col_grid                                 ! grid with time collection
            type(XML_str_type) :: top                                           ! topology
            type(XML_str_type) :: XYZ(3)                                        ! data items for geometry
            type(XML_str_type) :: geom                                          ! geometry
            integer :: loc_dim(3,2)                                             ! local dimensions (flux surfaces, field lines)
            integer :: n_r                                                      ! nr. of normal points
            
            ! initialize ierr
            ierr = 0
            
            ! only local masters
            if (rank.eq.0) then
                ! user output
                call writo('Drawing animation with HDF5')
                call lvl_ud(1)
                
                ! set up loc_dim and n_r
                loc_dim(:,1) = [size(X_1,1),size(X_1,2),1]
                loc_dim(:,2) = [size(X_2,1),size(X_2,2),1]
                n_r = size(X_1,3) - 1                                           ! should be same for all other X_i, Y_i and Z_i
                
                ! set up plot titles
                plot_title(1) = 'Magnetic Flux Surfaces'
                plot_title(2) = 'Magnetic Field Lines'
                
                ! open HDF5 file
                ierr = open_HDF5_file(file_info,'field_lines_in_flux_surfaces',&
                    &description=anim_name,ind_plot=.true.)
                CHCKERR('')
                
                ! create grid for time collection
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
                    ierr = print_HDF5_3D_data_item(XYZ(3),file_info,'Z_V_surf_'&
                        &//trim(i2str(id)),Z_1(:,:,id:id),loc_dim(:,1))
                    CHCKERR('')
                    
                    ! print geometry with X, Y and Z data item
                    call print_HDF5_geom(geom,2,XYZ,.true.)
                    
                    ! create a grid with the topology and the geometry
                    ierr = print_HDF5_grid(grids(1),plot_title(1),1,&
                        &grid_top=top,grid_geom=geom,reset=.true.)
                    CHCKERR('')
                    
                    ! B. end with magnetic field
                    
                    ! print topology
                    call print_HDF5_top(top,2,loc_dim(:,2))
                    
                    ! print data item for X
                    ierr = print_HDF5_3D_data_item(XYZ(1),file_info,'X_field_'&
                        &//trim(i2str(id)),X_2(:,:,id+1:id+1),loc_dim(:,2))
                    CHCKERR('')
                    
                    ! print data item for Y
                    ierr = print_HDF5_3D_data_item(XYZ(2),file_info,'Y_field_'&
                        &//trim(i2str(id)),Y_2(:,:,id+1:id+1),loc_dim(:,2))
                    CHCKERR('')
                    
                    ! print data item for Z
                    ierr = print_HDF5_3D_data_item(XYZ(3),file_info,'Z_field_'&
                        &//trim(i2str(id)),Z_2(:,:,id+1:id+1),loc_dim(:,2))
                    CHCKERR('')
                    
                    ! print geometry with X, Y and Z data item
                    call print_HDF5_geom(geom,2,XYZ,.true.)
                    
                    ! create a grid with the topology and the geometry
                    ierr = print_HDF5_grid(grids(2),plot_title(2),1,&
                        &grid_top=top,grid_geom=geom,reset=.true.)
                    CHCKERR('')
                    
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
                call add_HDF5_item(file_info,time_col_grid,reset=.true.)
                
                ! close HDF5 file
                ierr = close_HDF5_file(file_info)
                CHCKERR('')
                
                call lvl_ud(-1)
            end if
        end function plot_grid_real_HDF5
    end function plot_grid_real
    
    ! Extend a  grid angularly using  equidistant variables of  n_theta_plot and
    ! n_zeta_plot angular and  own loc_n_r points in  E coordinates. Optionally,
    ! the grid can  also be converted to  F coordinates if equilibrium  grid and
    ! eq variables are provided.
    integer function extend_grid_E(grid_in,grid_ext,grid_eq,eq) result(ierr)
        use num_vars, only: n_theta_plot, n_zeta_plot
        use grid_vars, only: create_grid
        
        character(*), parameter :: rout_name = 'extend_grid_E'
        
        ! input / output
        type(grid_type), intent(in) :: grid_in                                  ! grid to be extended
        type(grid_type), intent(inout) :: grid_ext                              ! extended grid
        type(grid_type), intent(in), optional :: grid_eq                        ! equilibrium grid
        type(eq_type), intent(in), optional :: eq                               ! equilibirum variables
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! tests
        if (present(grid_eq).neqv.present(eq)) then
            ierr = 1
            err_msg = 'When converting to Flux coordinates, also equilibrium &
                &grid and variables needed'
            CHCKERR(err_msg)
        end if
        
        ! creating  equilibrium  grid  for  the output  that  covers  the  whole
        ! geometry angularly in E coordinates
        ierr = create_grid(grid_ext,&
            &[n_theta_plot,n_zeta_plot,grid_in%n(3)],&
            &[grid_in%i_min,grid_in%i_max])
        CHCKERR('')
        grid_ext%r_E = grid_in%r_E
        grid_ext%loc_r_E = grid_in%loc_r_E
        ierr = calc_eqd_grid(grid_ext%theta_E,1*pi,3*pi,1)                      ! starting from pi gives nicer plots
        CHCKERR('')
        ierr = calc_eqd_grid(grid_ext%zeta_E,0*pi,2*pi,2)
        CHCKERR('')
        
        ! convert all E coordinates to F coordinates if requested
        if (present(grid_eq)) then
            grid_ext%r_F = grid_in%r_F
            ierr = coord_E2F(eq,grid_eq,&
                &grid_ext%loc_r_E,grid_ext%theta_E,grid_ext%zeta_E,&
                &grid_ext%loc_r_F,grid_ext%theta_F,grid_ext%zeta_F)
            CHCKERR('')
        end if
    end function extend_grid_E
    
    ! calculates loc_r_E (style 1) or loc_r_F (style 2)
    integer function calc_loc_r(grid_eq,eq,loc_r,style) result (ierr)
        use num_vars, only: use_pol_flux_E, use_pol_flux_F, eq_style
        use eq_vars, only: max_flux_p_F, max_flux_t_F
        
        character(*), parameter :: rout_name = 'calc_loc_r'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid
        type(eq_type), intent(in) :: eq                                         ! equilibrium variables
        real(dp), intent(inout), allocatable :: loc_r(:)                        ! loc_r
        integer, intent(in) :: style                                            ! whether to calculate in E or F
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! allocate loc_r if necessary
        if (allocated(loc_r)) then
            if (size(loc_r).ne.grid_eq%loc_n_r) then
                ierr = 1
                err_msg = 'loc_r needs to have the correct dimensions'
                CHCKERR(err_msg)
            end if
        else
            allocate(loc_r(grid_eq%loc_n_r))
        end if
        
        ! choose which style is being used:
        select case (style)
            case (1)                                                            ! E coords.
                ! choose which equilibrium style is being used:
                !   1:  VMEC
                !   2:  HELENA
                select case (eq_style)
                    case (1)                                                    ! VMEC
                        if (use_pol_flux_E) then                                ! normalized poloidal flux
                            loc_r = eq%flux_p_FD(:,0)/max_flux_p_F
                        else                                                    ! normalized toroidal flux
                            loc_r = eq%flux_t_FD(:,0)/max_flux_t_F
                        end if
                    case (2)                                                    ! HELENA
                        loc_r = eq%flux_p_FD(:,0)/(2*pi)                        ! poloidal flux / 2pi
                    case default
                        err_msg = 'No equilibrium style associated with '//&
                            &trim(i2str(eq_style))
                        ierr = 1
                        CHCKERR(err_msg)
                end select
            case (2)                                                            ! F coords.
                if (use_pol_flux_F) then                                        ! poloidal flux / 2pi
                    loc_r = eq%flux_p_FD(:,0)/(2*pi)
                else                                                            ! toroidal flux / 2pi
                    loc_r = eq%flux_t_FD(:,0)/(2*pi)
                end if
            case default
                err_msg = 'No style associated with '//trim(i2str(style))
                ierr = 1
                CHCKERR(err_msg)
        end select
    end function calc_loc_r
    
    ! Calculates volume integral on a 3D grid.
    ! Two angular and  one normal variable has  to be provided on a the grid. If
    ! the  i'th dimension  of the  grid  is equal  to  one, the  function to  be
    ! integrated is assumed not to vary in this dimension.
    ! Furthermore, if i is 1 or  2, the corresponding i'th (angular) variable is
    ! the only  variable that is  assumed to vary  in that dimension.  The other
    ! angular variable as well as the normal variable are assumed to be constant
    ! like the  function itself. However,  if i is 3,  an error is  displayed as
    ! this does not represent a physical situation.
    ! A common  case through which to  understand this is the  axisymmetric case
    ! where the first  angular variable theta varies in the  dimensions 1 and 3,
    ! the second angular  variable zeta varies only in the  dimension 2, and the
    ! normal variable only varies in the dimension 3.
    ! Alternatively,  there is  the case  of a  grid-aligned set  of coordinates
    ! theta and  alpha, where the  first dimension corresponds to  the direction
    ! along the magnetic field line, the  second to the geodesical direction and
    ! the thirs to the normal direction. If the calculations for different field
    ! lines are  decoupled, the variation in  the second dimension is  not taken
    ! into account and no integration happens along it.
    ! Internally, the angular  variables and the normal variable  are related to
    ! the coordinates (x,y,z) that correspond to the three dimensions. They thus
    ! form a computational orthogonal grid to which the original coordinates are
    ! related through the transformation of Jacobians:
    !   Jxyz = J r_F_z (theta_x zeta_y - theta_y zeta_x) ,
    ! so that the integral becomes
    !   sum_xyz f(x,y,z) J(x,y,z) r_F_z (theta_x zeta_y - theta_y zeta_x) dxdydz
    ! where dx, dy and dz are all trivially equal to 1.
    ! The integrand has to be evaluated at the intermediate positions inside the
    ! cells. This is  done by taking the  average of the 2^3=8 points  for fJ as
    ! well as the transformation of the Jacobian.
    integer function calc_int_vol(ang_1,ang_2,norm,J,f,f_int) result(ierr)
#if ldebug
        use num_vars, only: rank
#endif
        
        character(*), parameter :: rout_name = 'calc_int_vol'
        
        ! input / output
        real(dp), intent(in) :: ang_1(:,:,:), ang_2(:,:,:), norm(:)             ! coordinate variables
        real(dp), intent(in) :: J(:,:,:)                                        ! Jacobian
        complex(dp), intent(in) :: f(:,:,:,:)                                   ! input f(n_par,n_geo,n_r,size_X^2)
        complex(dp), intent(inout) :: f_int(:)                                  ! output integrated f
        
        ! local variables
        integer :: id, jd, kd, ld                                               ! counters
        integer :: nn_mod                                                       ! number of indices for V and V_int
        integer :: k_min                                                        ! minimum k
        integer :: dims(3)                                                      ! real dimensions
        real(dp), allocatable :: transf_J(:,:,:,:)                              ! comps. of transf. between J and Jxyz: theta_x, zeta_y, theta_y, zeta_x and r_F_z
        real(dp), allocatable :: transf_J_tot(:,:,:)                            ! transf. between J and Jxyz: (theta_x*zeta_y-theta_y*zeta_x)*r_F_z
        complex(dp), allocatable :: Jf(:,:,:,:)                                 ! J*f
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer :: i_a(3), i_b(3)                                               ! limits of building blocks of intermediary variables
        logical :: dim_1(2)                                                     ! whether the dimension is equal to one
#if ldebug
        character(len=max_str_ln), allocatable :: var_names(:,:)                ! names of variables
#endif
        
        ! initialize ierr
        ierr = 0
        
        ! set nn_mod
        nn_mod = size(f,4)
        
        ! tests
        if (size(f_int).ne.nn_mod) then
            ierr = 1
            err_msg = 'f and f_int need to have the same storage convention'
            CHCKERR(err_msg)
        end if
        if (size(ang_1).ne.size(ang_2) .or. size(ang_1,3).ne.size(norm)) then
            ierr = 1
            err_msg = 'The angular variables have to be compatible'
            CHCKERR(err_msg)
        end if
        if (size(norm).eq.1) then
            ierr = 1
            err_msg = 'The normal grid size has to be greater than one'
            CHCKERR(err_msg)
        end if
        
        ! set up dims
        dims = shape(ang_1)
        
        ! set up dim_1
        dim_1 = .false.
        if (size(ang_1,1).eq.1) dim_1(1) = .true.
        if (size(ang_2,2).eq.1) dim_1(2) = .true.
        
        ! set up Jf and transf_J
        allocate(Jf(max(dims(1)-1,1),max(dims(2)-1,1),dims(3)-1,nn_mod))
        allocate(transf_J(max(dims(1)-1,1),max(dims(2)-1,1),dims(3)-1,5))
        allocate(transf_J_tot(max(dims(1)-1,1),max(dims(2)-1,1),dims(3)-1))
        
        ! set up k_min
        k_min = 1
        
        ! get intermediate Jf and components of transf_J:
        !   1: ang_1_x
        !   2: ang_2_y
        !   3: ang_1_y
        !   4: ang_2_x
        !   5: r_F_z
        ! Note: the are always 8 terms in the summation, hence the factor 0.125
        Jf = 0._dp
        transf_J = 0._dp
        do kd = 1,2
            if (kd.eq.1) then
                i_a(3) = 1
                i_b(3) = dims(3)-1
            else
                i_a(3) = 2
                i_b(3) = dims(3)
            end if
            do jd = 1,2
                if (dim_1(2)) then
                        i_a(2) = 1
                        i_b(2) = dims(2)
                else
                    if (jd.eq.1) then
                        i_a(2) = 1
                        i_b(2) = dims(2)-1
                    else
                        i_a(2) = 2
                        i_b(2) = dims(2)
                    end if
                end if
                do id = 1,2
                    if (dim_1(1)) then
                            i_a(1) = 1
                            i_b(1) = dims(1)
                    else
                        if (id.eq.1) then
                            i_a(1) = 1
                            i_b(1) = dims(1)-1
                        else
                            i_a(1) = 2
                            i_b(1) = dims(1)
                        end if
                    end if
                    do ld = 1,nn_mod
                        ! Jf
                        Jf(:,:,:,ld) = Jf(:,:,:,ld) + 0.125_dp * &
                            &J(i_a(1):i_b(1),i_a(2):i_b(2),i_a(3):i_b(3)) * &
                            &f(i_a(1):i_b(1),i_a(2):i_b(2),i_a(3):i_b(3),ld)
                    end do
                    if (dim_1(1)) then
                        ! transf_J(1): ang_1_x
                        transf_J(:,:,:,1) = 1._dp
                        ! transf_J(4): ang_2_x
                        transf_J(:,:,:,4) = 0._dp
                    else
                        do ld = 1,dims(1)-1
                            ! transf_J(1): ang_1_x
                            transf_J(ld,:,:,1) = transf_J(ld,:,:,1) + &
                                &0.125_dp * ( &
                                &ang_1(ld+1,i_a(2):i_b(2),i_a(3):i_b(3)) - &
                                &ang_1(ld,i_a(2):i_b(2),i_a(3):i_b(3)) )
                            ! transf_J(4): ang_2_x
                            transf_J(ld,:,:,4) = transf_J(ld,:,:,4) + &
                                &0.125_dp * ( &
                                &ang_2(ld+1,i_a(2):i_b(2),i_a(3):i_b(3)) - &
                                &ang_2(ld,i_a(2):i_b(2),i_a(3):i_b(3)) )
                        end do
                    end if
                    if (dim_1(2)) then
                        ! transf_J(2): ang_2_y
                        transf_J(:,:,:,2) = 1._dp
                        ! transf_J(3): ang_1_y
                        transf_J(:,:,:,3) = 0._dp
                    else
                        do ld = 1,dims(2)-1
                            ! transf_J(2): ang_2_y
                            transf_J(:,ld,:,2) = transf_J(:,ld,:,2) + &
                                &0.125_dp * ( &
                                &ang_2(i_a(1):i_b(1),ld+1,i_a(3):i_b(3)) - &
                                &ang_2(i_a(1):i_b(1),ld,i_a(3):i_b(3)) )
                            ! transf_J(3): ang_1_y
                            transf_J(:,ld,:,3) = transf_J(:,ld,:,3) + &
                                &0.125_dp * ( &
                                &ang_1(i_a(1):i_b(1),ld+1,i_a(3):i_b(3)) - &
                                &ang_1(i_a(1):i_b(1),ld,i_a(3):i_b(3)) )
                        end do
                    end if
                    do ld = 1,dims(3)-1
                        ! transf_J(5): r_F_z
                        transf_J(:,:,ld,5) = transf_J(:,:,ld,5) + 0.125_dp * ( &
                            &norm(ld+1) - norm(ld) )
                    end do
                end do
            end do
        end do
        
        ! set up total transf_J
        transf_J_tot = (transf_J(:,:,:,1)*transf_J(:,:,:,2)-&
            &transf_J(:,:,:,3)*transf_J(:,:,:,4))*transf_J(:,:,:,5)
        
        ! integrate term  over ang_par_F for all equilibrium  grid points, using
        ! dx = dy = dz =1
        do ld = 1,nn_mod
            f_int(ld) = sum(Jf(:,:,:,ld)*transf_J_tot)
        end do
        
#if ldebug
        if (debug_calc_int_vol) then
            call writo('Testing whether volume integral is calculated &
                &correctly')
            call lvl_ud(1)
            
            allocate(var_names(nn_mod,2))
            var_names(:,1) = 'real part of integrand J f'
            var_names(:,2) = 'imaginary part of integrand J f'
            do ld = 1,nn_mod
                var_names(ld,1) = trim(var_names(ld,1))//' '//trim(i2str(ld))
                var_names(ld,2) = trim(var_names(ld,2))//' '//trim(i2str(ld))
            end do
            
            call plot_HDF5(var_names(:,1),'TEST_RE_Jf_'//trim(i2str(rank)),&
                &realpart(Jf))
            call plot_HDF5(var_names(:,2),'TEST_IM_Jf_'//trim(i2str(rank)),&
                &imagpart(Jf))
            call plot_HDF5('transformation of Jacobians',&
                &'TEST_transf_J_'//trim(i2str(rank)),transf_J_tot)
            
            call writo('Resulting integral: ')
            call lvl_ud(1)
            do ld = 1,nn_mod
                call writo('integral '//trim(i2str(ld))//': '//&
                    &trim(c2str(f_int(ld))))
            end do
            call lvl_ud(-1)
            
            call lvl_ud(-1)
        end if
#endif
        
        ! deallocate local variables
        deallocate(Jf,transf_J,transf_J_tot)
    end function calc_int_vol
    
    ! calculate interpolation  factors for  normal interpolation in  grid_out of
    ! quantities defined on grid_in.
    ! The output  loc_r consists of  an array of  real values where  the floored
    ! integer of each  value indicates the base index of  the interpolated value
    ! in the output  grid and the modulus  is the the fraction  towards the next
    ! integer.
    ! By default, the variables in the Flux coord. system are used, but this can
    ! be changed optionally with the flag "use_E".
    integer function get_norm_interp_data(grid_in,grid_out,loc_r,use_E) &
        &result(ierr)
        use utilities, only: con2dis
        
        character(*), parameter :: rout_name = 'get_norm_interp_data'
        
        ! input / output
        type(grid_type), intent(in) :: grid_in, grid_out                        ! input and output grid
        real(dp), allocatable, intent(inout) :: loc_r(:)                        ! interpolation index
        logical, intent(in), optional :: use_E                                  ! whether E is used instead of F
        
        ! local variables
        integer :: kd                                                           ! counter
        integer :: loc_n_r                                                      ! nr. of normal points in output grid
        logical :: use_E_loc                                                    ! local version of use_E
        real(dp), pointer :: loc_r_in(:) => null()                              ! loc_r_F or loc_r_E of grids
        real(dp), pointer :: loc_r_out(:) => null()                             ! loc_r_F or loc_r_E of grids
#if ldebug
        real(dp), allocatable :: loc_r_ALT(:)                                   ! alternative calculation for loc_r
#endif
        
        ! initialize ierr
        ierr = 0
        
        ! set up local use_E
        use_E_loc = .false.
        if (present(use_E)) use_E_loc = use_E
        
        ! set loc_n_r
        loc_n_r = size(grid_out%loc_r_F)
        
        ! allocate loc_r
        allocate(loc_r(loc_n_r))
        
        ! point loc_r_in and loc_r_out
        if (use_E_loc) then
            loc_r_in => grid_in%loc_r_E
            loc_r_out => grid_out%loc_r_E
        else
            loc_r_in => grid_in%loc_r_F
            loc_r_out => grid_out%loc_r_F
        end if
        
#if ldebug
        if (debug_get_norm_interp_data) then
            call print_GP_2D('loc_r_in','',loc_r_in)
            call print_GP_2D('loc_r_out','',loc_r_out)
        end if
#endif
        
        ! get the table  indices between which to interpolate  by converting the
        ! continuous equilibrium points to the discretized values
        ! loop over all normal points
        do kd = 1,loc_n_r
            ierr = con2dis(loc_r_out(kd),loc_r(kd),loc_r_in)
            CHCKERR('')
        end do
        
#if ldebug
        if (debug_get_norm_interp_data) then
            ! now use the other variables to see whether the result is the same
            ! point loc_r_in and loc_r_out
            if (.not.use_E_loc) then
                loc_r_in => grid_in%loc_r_E
                loc_r_out => grid_out%loc_r_E
            else
                loc_r_in => grid_in%loc_r_F
                loc_r_out => grid_out%loc_r_F
            end if
            ! allocate alternative loc_r
            allocate(loc_r_ALT(loc_n_r))
            ! loop over all normal points
            do kd = 1,loc_n_r
                ierr = con2dis(loc_r_out(kd),loc_r_ALT(kd),loc_r_in)
                CHCKERR('')
            end do
            ! plot output
            !call print_GP_2D('loc_r','',loc_r)
            !call print_GP_2D('ALTERNATIVE loc_r','',loc_r_ALT)
            call writo('Maximum difference: '//&
                &trim(r2str(maxval(abs(loc_r-loc_r_ALT)))))
        end if
#endif
        
        ! clean up
        nullify(loc_r_in,loc_r_out)
    end function get_norm_interp_data
    
    ! Trim a grid, removing any overlap between the different regions.
    ! By default, the routine assumes a  symmetric ghost region and cuts as many
    ! grid points from the end of the  previous process as from the beginning of
    ! the next process, but if the number of overlapping grid points is odd, the
    ! previous process looses one more point.
    ! Optionally, the trimmed indices in the normal direction can be provided.
    integer function trim_grid(grid_in,grid_out,norm_id) result(ierr)
        use num_vars, only: n_procs, rank
        use MPI_utilities, only: get_ser_var
        use grid_vars, only: create_grid
#if ldebug
        use MPI_utilities, only: wait_MPI
#endif
        
        character(*), parameter :: rout_name = 'trim_grid'
        
        ! input / output
        type(grid_type), intent(in) :: grid_in                                  ! input grid
        type(grid_type), intent(inout) :: grid_out                              ! trimmed grid
        integer, intent(inout), optional :: norm_id(2)                          ! normal indices corresponding to trimmed part
        
        ! local variables
        integer, allocatable :: tot_i_min(:)                                    ! i_min of grid of all processes
        integer, allocatable :: tot_i_max(:)                                    ! i_max of grid of all processes
        integer :: i_lim_out(2)                                                 ! i_lim of output grid
        integer :: n_out(3)                                                     ! n of output grid
        integer :: norm_shift                                                   ! shift of normal indices
        
        ! initialize ierr
        ierr = 0
        
        ! detect whether grid divided
        if (grid_in%divided) then
            ! get min_i's of the grid_in
            ierr = get_ser_var([grid_in%i_min],tot_i_min,scatter=.true.)
            CHCKERR('')
            
            ! get max_i's of the grid_in
            ierr = get_ser_var([grid_in%i_max],tot_i_max,scatter=.true.)
            CHCKERR('')
            
            ! set  i_lim of trimmed output  grid (not yet shifted  by first proc
            ! min)
            if (rank.gt.0) then
                i_lim_out(1) = max(tot_i_min(1),tot_i_min(rank+1)+&
                    &floor((tot_i_max(rank)-tot_i_min(rank+1)+1._dp)/2))
            else
                i_lim_out(1) = tot_i_min(1)
            end if
            if (rank.lt.n_procs-1) then
                i_lim_out(2) = min(tot_i_max(n_procs),tot_i_max(rank+1)-&
                    &ceiling((tot_i_max(rank+1)-tot_i_min(rank+2)+1._dp)/2))
            else
                i_lim_out(2) = tot_i_max(n_procs)
            end if
            
            ! normal shift between grids
            norm_shift = i_lim_out(1) - grid_in%i_min
            
            ! get  min_i's and max_i's  of the grid_out,  not shifted by  min of
            ! first process
            ierr = get_ser_var([i_lim_out(1)],tot_i_min,scatter=.true.)
            CHCKERR('')
            ierr = get_ser_var([i_lim_out(2)],tot_i_max,scatter=.true.)
            CHCKERR('')
            
            ! set n of output grid
            n_out(1) = grid_in%n(1)
            n_out(2) = grid_in%n(2)
            n_out(3) = sum(tot_i_max-tot_i_min+1)
            
            ! create new grid
            ierr = create_grid(grid_out,n_out,i_lim_out-tot_i_min(1)+1)         ! limits shifted by min of first process
            CHCKERR('')
            
            ! recycle  i_lim_out  for  shifted  array  indices, set  norm_id  if
            ! requested
            i_lim_out = i_lim_out - i_lim_out(1) + 1 + norm_shift
            if (present(norm_id)) norm_id = i_lim_out
        else
            ! set n of output grid
            n_out = grid_in%n
            
            ! set unshifted i_lim_out and norm_id if requested
            i_lim_out = [grid_in%i_min,grid_in%i_max]
            if (present(norm_id)) norm_id = i_lim_out
            
            ! shift i_lim_out by min.
            i_lim_out = i_lim_out-i_lim_out(1)+1
            
            ! create new grid
            ierr = create_grid(grid_out,n_out,i_lim_out)                        ! grid not divided
            CHCKERR('')
        end if
        
        ! copy local arrays
        if (grid_in%n(1).ne.0 .and. grid_in%n(2).ne.0) then                     ! only if 3D grid
            grid_out%theta_E = grid_in%theta_E(:,:,i_lim_out(1):i_lim_out(2))
            grid_out%zeta_E = grid_in%zeta_E(:,:,i_lim_out(1):i_lim_out(2))
            grid_out%theta_F = grid_in%theta_F(:,:,i_lim_out(1):i_lim_out(2))
            grid_out%zeta_F = grid_in%zeta_F(:,:,i_lim_out(1):i_lim_out(2))
        end if
        if (grid_in%divided) then                                               ! but if input grid divided, loc_r gets priority
            grid_out%loc_r_E = grid_in%loc_r_E(i_lim_out(1):i_lim_out(2))
            grid_out%loc_r_F = grid_in%loc_r_F(i_lim_out(1):i_lim_out(2))
        end if
        
        ! if divided, set total arrays
        if (grid_in%divided) then
            grid_out%r_E = grid_in%r_E(tot_i_min(1):tot_i_max(n_procs))
            grid_out%r_F = grid_in%r_F(tot_i_min(1):tot_i_max(n_procs))
        else
            grid_out%r_E = grid_in%r_E
            grid_out%r_F = grid_in%r_F
        end if
    end function trim_grid
    
    ! Untrims a trimmed  grid by introducing an assymetric ghost  regions at the
    ! right. The width of the ghost region has to be provided.
    ! Note: The ghosted grid should be deallocated (with dealloc_grid).
    ! Note: The input grid HAS to be trimmed!
    integer function untrim_grid(grid_in,grid_out,size_ghost) result(ierr)
        use num_vars, only: n_procs, rank
        use MPI_utilities, only: get_ghost_arr, get_ser_var
        use grid_vars, only: create_grid
        
        character(*), parameter :: rout_name = 'untrim_grid'
        
        ! input / output
        type(grid_type), intent(in) :: grid_in                                  ! input grid
        type(grid_type), intent(inout) :: grid_out                              ! ghosted grid
        integer, intent(in) :: size_ghost                                       ! width of ghost region
        
        ! local variables
        integer :: i_lim_in(2)                                                  ! limits of input grid
        integer :: i_lim_out(2)                                                 ! limits of ghosted grid
        integer, allocatable :: tot_loc_n_r(:)                                  ! loc_n_r of all processes
        integer :: size_ghost_loc                                               ! local size_ghost
        
        ! initialize ierr
        ierr = 0
        
        ! get array sizes of all processes
        ierr = get_ser_var([grid_in%loc_n_r],tot_loc_n_r,scatter=.true.)
        CHCKERR('')
        
        ! set local size_ghost
        size_ghost_loc = min(size_ghost,minval(tot_loc_n_r)-1)
        
        ! set i_lim_in and i_lim_out
        i_lim_in = [grid_in%i_min,grid_in%i_max]
        i_lim_out = i_lim_in
        if (rank.lt.n_procs-1) i_lim_out(2) = i_lim_out(2)+size_ghost_loc
        
        ! create grid
        ierr = create_grid(grid_out,grid_in%n,i_lim_out)
        CHCKERR('')
        
        ! set ghosted variables
        if (grid_in%n(1).ne.0 .and. grid_in%n(2).ne.0) then                     ! only if 3D grid
            grid_out%theta_E(:,:,1:grid_in%loc_n_r) = grid_in%theta_E
            grid_out%zeta_E(:,:,1:grid_in%loc_n_r) = grid_in%zeta_E
            grid_out%theta_F(:,:,1:grid_in%loc_n_r) = grid_in%theta_F
            grid_out%zeta_F(:,:,1:grid_in%loc_n_r) = grid_in%zeta_F
            ierr = get_ghost_arr(grid_out%theta_E,size_ghost_loc)
            CHCKERR('')
            ierr = get_ghost_arr(grid_out%zeta_E,size_ghost_loc)
            CHCKERR('')
            ierr = get_ghost_arr(grid_out%theta_F,size_ghost_loc)
            CHCKERR('')
            ierr = get_ghost_arr(grid_out%zeta_F,size_ghost_loc)
            CHCKERR('')
        end if
        if (grid_in%divided) then                                               ! but if input grid divided, loc_r gets priority
            grid_out%loc_r_E = grid_in%r_E(i_lim_out(1):i_lim_out(2))
            grid_out%loc_r_F = grid_in%r_F(i_lim_out(1):i_lim_out(2))
        end if
        grid_out%r_E = grid_in%r_E
        grid_out%r_F = grid_in%r_F
    end function untrim_grid
end module grid_ops
