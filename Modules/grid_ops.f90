!------------------------------------------------------------------------------!
!   Operations that have to do with the grids and different coordinate systems !
!------------------------------------------------------------------------------!
module grid_ops
#include <PB3D_macros.h>
    use str_ops
    use output_ops
    use messages
    use num_vars, only: dp, pi, max_str_ln
    use grid_vars, only: grid_type, disc_type
    use eq_vars, only: eq_1_type

    implicit none
    private
    public calc_norm_range, calc_ang_grid_eq_B, magn_grid_plot, setup_grid_eq, &
        &setup_grid_eq_B, setup_grid_X, setup_grid_sol, print_output_grid
#if ldebug
    public debug_calc_ang_grid_eq_B
#endif
    
    ! global variables
#if ldebug
    logical :: debug_calc_ang_grid_eq_B = .false.                               ! plot debug information for calc_ang_grid_eq_B
#endif

contains
    ! Calculates  normal range for the  input grid, the equilibrium  grid and/or
    ! the solution grid.
    ! For PB3D:
    !   - if  'in_limits' is provided, the  limits are found by  calculating the
    !   tightest equilibrium grid points that  encompass the solution range. The
    !   global max_flux variables are also set.
    !   - if  'eq_limits'  is provided,  the limits  are found  by dividing  the
    !   equilibrium range over the available processes. 'r_F_eq' is ignored.
    !   - if 'X_limits'  is provided, the limits  are set equal to  the case for
    !   'sol_limits'.
    !   - if 'sol_limits'  is provided, the solution range is  divided under the
    !   processes and 'r_F_eq' is filled. Note that  it has to be of the correct
    !   size for this.
    ! For POST:
    !   -  The  solution  range  is calculated  and  the  tightest  encompassing
    !   equilibrium range is found. Also 'r_F_eq' and 'r_F_sol' are provided and
    !   necessary.
    ! Note: The input  limits strip the input variables from  the ranges that do
    ! not  matter to  the  solution  range requested  but  does  not divide  the
    ! resulting range  under processes.  This is  the information  of eq_limits,
    ! which implies  that in_limits have  to be  correctly applied to  the input
    ! variables for eq_limits to behave correctly.
    integer function calc_norm_range(in_limits,eq_limits,X_limits,sol_limits,&
        &r_F_eq,r_F_X,r_F_sol) result(ierr)
        use num_vars, only: prog_style
        
        character(*), parameter :: rout_name = 'calc_norm_range'
        
        ! input / output
        integer, intent(inout), optional :: in_limits(2)                        ! min. and max. index of in grid
        integer, intent(inout), optional :: eq_limits(2)                        ! min. and max. index of eq grid for this process
        integer, intent(inout), optional :: X_limits(2)                         ! min. and max. index of X grid for this process
        integer, intent(inout), optional :: sol_limits(2)                       ! min. and max. index of sol grid for this process
        real(dp), intent(inout), optional :: r_F_eq(:)                          ! equilibrium r_F
        real(dp), intent(inout), optional :: r_F_X(:)                           ! perturbation r_F
        real(dp), intent(inout), optional :: r_F_sol(:)                         ! solution r_F
        
        ! initialize ierr
        ierr = 0
        
        ! select depending on program style
        select case (prog_style)
            case(1)                                                             ! PB3D
                if (present(in_limits)) then
                    ierr = calc_norm_range_PB3D_in(in_limits)
                    CHCKERR('')
                end if
                if (present(eq_limits)) then
                    call calc_norm_range_PB3D_eq(eq_limits)
                end if
                if (present(X_limits).and.present(r_F_X)) then
                    ierr = calc_norm_range_PB3D_X(X_limits,r_F_X)
                    CHCKERR('')
                end if
                if (present(sol_limits).and.present(r_F_sol)) then
                    ierr = calc_norm_range_PB3D_sol(sol_limits,r_F_sol)
                    CHCKERR('')
                end if
            case(2)                                                             ! PB3D post-processing
                if (.not.present(eq_limits) .or. .not.present(X_limits) .or. &
                    &.not.present(sol_limits) .or. .not.present(r_F_eq) .or. &
                    &.not. present(r_F_sol)) then
                    ierr = 1
                    CHCKERR('Incorrect variables provided.')
                end if
                ierr = calc_norm_range_POST(eq_limits,X_limits,sol_limits,&
                    &r_F_eq,r_F_sol)
                CHCKERR('')
        end select
    contains
        ! The normal  range is calculated by  finding the tightest range  of the
        ! input variables encompassing the  entire solution range, including the
        ! tolerance "tol_norm". Additionally, the  global max_flux variables are
        ! set as well..
        ! Note: this is the global, undivided range. The division information is
        ! calculated in 'eq_limits'. Furthermore, the input variables have to be
        ! tabulated on the full grid provided by the equilibrium code.
        integer function calc_norm_range_PB3D_in(in_limits) &
            &result(ierr)                                                       ! PB3D version for equilibrium grid
            use num_vars, only: use_pol_flux_E, use_pol_flux_F, eq_style, &
                &tol_norm, norm_disc_prec_eq
            use num_utilities, only: con2dis, dis2con, calc_int, round_with_tol
            use grid_utilities, only: setup_deriv_data, setup_interp_data, &
                &apply_disc
            use VMEC, only: flux_p_V, flux_t_V
            use HELENA_vars, only: flux_p_H, qs_H
            use X_vars, only: min_r_sol, max_r_sol
            use eq_vars, only: max_flux_E, max_flux_F
            use grid_vars, only: n_r_in
            
            character(*), parameter :: rout_name = 'calc_norm_range_PB3D_in'
            
            ! input / output
            integer, intent(inout) :: in_limits(2)                              ! total min. and max. index of eq. grid for this process
            
            ! local variables
            real(dp), allocatable :: flux_F(:), flux_E(:)                       ! either pol. or tor. flux in F and E
            real(dp), allocatable :: Dflux_p_H(:)                               ! normal derivative of flux_p_H
            type(disc_type) :: norm_deriv_data                                  ! data for normal derivatives
            real(dp) :: tot_min_r_in_F_con                                      ! tot_min_r_in in continuous F coords.
            real(dp) :: tot_max_r_in_F_con                                      ! tot_max_r_in in continuous F coords.
            real(dp) :: tot_min_r_in_E_con(1)                                   ! tot_min_r_in in continuous E coords.
            real(dp) :: tot_max_r_in_E_con(1)                                   ! tot_max_r_in in continuous E coords.
            real(dp) :: tot_min_r_in_E_dis                                      ! tot_min_r_in in discrete E coords., unrounded
            real(dp) :: tot_max_r_in_E_dis                                      ! tot_max_r_in in discrete E coords., unrounded
            type(disc_type) :: norm_interp_data                                 ! data for normal interpolation
            
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
                        flux_F = flux_p_V
                    else
                        flux_F = - flux_t_V                                     ! conversion VMEC LH -> RH coord. system
                    end if
                    ! set up E flux
                    if (use_pol_flux_E) then
                        flux_F = flux_p_V
                    else
                        flux_E = flux_t_V
                    end if
                case (2)                                                        ! HELENA
                    ! calculate normal derivative of flux_p_H
                    allocate(Dflux_p_H(n_r_in))
                    ierr = setup_deriv_data(flux_p_H/(2*pi),norm_deriv_data,&
                        &1,norm_disc_prec_eq)
                    CHCKERR('')
                    ierr = apply_disc(flux_p_H,norm_deriv_data,Dflux_p_H)
                    CHCKERR('')
                    call norm_deriv_data%dealloc()
                    ! set up F flux
                    if (use_pol_flux_F) then
                        flux_F = flux_p_H
                    else
                        ierr = calc_int(qs_H*Dflux_p_H,flux_p_H/(2*pi),flux_F)
                        CHCKERR('')
                    end if
                    ! set up E flux
                    if (use_pol_flux_E) then
                        flux_E = flux_p_H
                    else
                        ierr = calc_int(qs_H*Dflux_p_H,flux_p_H/(2*pi),flux_E)
                        CHCKERR('')
                    end if
            end select
            
            ! set max_flux
            max_flux_E = flux_E(n_r_in)
            max_flux_F = flux_F(n_r_in)
            
            ! normalize flux_F and flux_E to (0..1) using max_flux
            flux_E = flux_E/max_flux_E
            flux_F = flux_F/max_flux_F
            
            ! get lower and upper bound of total solution range
            ! 1. include tolerance
            tot_min_r_in_F_con = max(min_r_sol-tol_norm,0._dp)
            tot_max_r_in_F_con = min(max_r_sol+tol_norm,1._dp)
            ! 2. round with tolerance
            ierr = round_with_tol(tot_min_r_in_F_con,0._dp,1._dp)
            CHCKERR('')
            ierr = round_with_tol(tot_max_r_in_F_con,0._dp,1._dp)
            CHCKERR('')
            ! 3. continuous E coords
            ierr = setup_interp_data(flux_F,[tot_min_r_in_F_con],&
                &norm_interp_data,norm_disc_prec_eq)
            CHCKERR('')
            ierr = apply_disc(flux_E,norm_interp_data,tot_min_r_in_E_con)
            CHCKERR('')
            ierr = setup_interp_data(flux_F,[tot_max_r_in_F_con],&
                &norm_interp_data,norm_disc_prec_eq)
            CHCKERR('')
            ierr = apply_disc(flux_E,norm_interp_data,tot_max_r_in_E_con)
            CHCKERR('')
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
            
            ! clean up
            call norm_interp_data%dealloc()
        end function calc_norm_range_PB3D_in
        
        ! The normal range is calculated  by dividing the full equilibrium range
        ! passed from the input phase between the processes.
        subroutine calc_norm_range_PB3D_eq(eq_limits)                           ! PB3D version for equilibrium grid
            use num_vars, only: n_procs, rank, norm_disc_prec_eq
            use grid_vars, only: n_r_eq
            
            ! input / output
            integer, intent(inout) :: eq_limits(2)                              ! min. and max. index of eq. grid for this process
            
            ! set local equilibrium limits
            eq_limits(1) = nint(1 + 1._dp*rank*(n_r_eq-1)/n_procs)
            eq_limits(2) = nint(1 + (1._dp+rank)*(n_r_eq-1)/n_procs)
            
            ! increase lower limits for processes that are not first
            if (rank.gt.0) eq_limits(1) = eq_limits(1)+1
            
            ! ghost regions of width 2*norm_disc_prec_eq
            eq_limits(1) = max(eq_limits(1)-2*norm_disc_prec_eq,1)
            eq_limits(2) = min(eq_limits(2)+2*norm_disc_prec_eq,n_r_eq)
        end subroutine calc_norm_range_PB3D_eq
        
        ! The normal range is determined by duplicating the normal range for the
        ! solution.  A divided  grid is  not necessary  as the  division in  the
        ! perturbation phase  is in the mode  numbers rather than in  the normal
        ! range.
        integer function calc_norm_range_PB3D_X(X_limits,r_F_X) &
            &result(ierr)                                                       ! PB3D version for perturbation grid
            character(*), parameter :: rout_name = 'calc_norm_range_PB3D_X'
            
            ! input / output
            integer, intent(inout) :: X_limits(2)                               ! min. and max. index of perturbation grid for this process
            real(dp), intent(inout) :: r_F_X(:)                                 ! perturbation r_F
            
            ! call the version for solution range
            ierr = calc_norm_range_PB3D_sol(X_limits,r_F_X)
            CHCKERR('')
            
            ! undivided grid
            X_limits = [1,size(r_F_X)]
        end function calc_norm_range_PB3D_X
        
        ! The normal range is determined  by simply dividing the solution range,
        ! no ghost range required.
        integer function calc_norm_range_PB3D_sol(sol_limits,r_F_sol) &
            &result(ierr)                                                       ! PB3D version for solution grid
            use num_vars, only: n_procs, rank
            use eq_vars, only: max_flux_F
            use X_vars, only: min_r_sol, max_r_sol
            use num_utilities, only: round_with_tol
            
            character(*), parameter :: rout_name = 'calc_norm_range_PB3D_sol'
            
            ! input / output
            integer, intent(inout) :: sol_limits(2)                             ! min. and max. index of sol grid for this process
            real(dp), intent(inout) :: r_F_sol(:)                               ! solution r_F
            
            ! local variables
            integer :: n_r_sol                                                  ! total nr. of normal points in solution grid
            integer :: id                                                       ! counter
            integer, allocatable :: loc_n_r_sol(:)                              ! local nr. of normal points in solution grid
            
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
            
            ! calculate r_F_sol in range from min_r_sol to max_r_sol
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
        integer function calc_norm_range_POST(eq_limits,X_limits,sol_limits,&
            &r_F_eq,r_F_sol) result(ierr)                                       ! POST version
            use num_vars, only: n_procs, rank, norm_disc_prec_sol
            
            character(*), parameter :: rout_name = 'calc_norm_range_POST'
            
            ! input / output
            integer, intent(inout) :: eq_limits(2), X_limits(2), sol_limits(2)  ! min. and max. index of eq, X and sol grid for this process
            real(dp), intent(in) :: r_F_eq(:), r_F_sol(:)                       ! eq and sol r_F
            
            ! local variables
            integer :: n_r_eq, n_r_sol                                          ! total nr. of normal points in eq and solution grid
            integer, allocatable :: loc_n_r_sol(:)                              ! local nr. of normal points in solution grid
            integer :: id                                                       ! counter
            real(dp) :: min_sol, max_sol                                        ! min. and max. of r_F_sol in range of this process
            real(dp), parameter :: tol = 1.E-6                                  ! tolerance for grids
            character(len=max_str_ln) :: err_msg                                ! error message
            
            ! initialize ierr
            ierr = 0
            
            ! initialize n_r_eq and n_r_sol
            n_r_eq = size(r_F_eq)
            n_r_sol = size(r_F_sol)
            allocate(loc_n_r_sol(n_procs))
            
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
            eq_limits = [0,n_r_eq+1]                                            ! initialize out of range
            if (r_F_eq(1).lt.r_F_eq(n_r_eq)) then                               ! ascending r_F_eq
                do id = 1,n_r_eq
                    if (r_F_eq(id).le.min_sol+tol) eq_limits(1) = id            ! move lower limit up
                    if (r_F_eq(n_r_eq-id+1).ge.max_sol-tol) &
                        &eq_limits(2) = n_r_eq-id+1                             ! move upper limit down
                end do
            else                                                                ! descending r_F_eq
                do id = 1,n_r_eq
                    if (r_F_eq(id).ge.max_sol+tol) eq_limits(1) = id            ! move lower limit up
                    if (r_F_eq(n_r_eq-id+1).le.min_sol-tol) &
                        &eq_limits(2) = n_r_eq-id+1                             ! move upper limit down
                end do
            end if
            
            ! check if valid limits found
            if (eq_limits(1).lt.1 .or. eq_limits(2).gt.n_r_eq) then
                ierr = 1
                err_msg = 'Solution range not contained in equilibrium range'
                CHCKERR(err_msg)
            end if
            
            ! copy solution range to perturbation range
            X_limits = sol_limits
        end function calc_norm_range_POST
    end function calc_norm_range

    ! Sets up the general equilibrium grid, in which the following variables are
    ! calculated:
    !   - equilibrium variables (eq)
    !   - perturbation variables (X)
    ! For the implementation of the equilibrium grid the normal part of the grid
    ! is always  given by the  output of the  equilibrium code, but  the angular
    ! part depends on the style:
    !   -  VMEC:  The  output  is   analytic  in  the  angular  coordinates,  so
    !   field-aligned coordinates  are used.  Later, in  calc_ang_grid_eq_B, the
    !   angular variables are calculated.
    !   - HELENA: The output is given on a poloidal grid (no toroidal dependency
    !   due to axisymmetry), which is used.
    ! Note: by  setting the flag "only_half_grid",  only the even points  of the
    ! parallel grid are set up. (Only for VMEC)
    integer function setup_grid_eq(grid_eq,eq_limits,only_half_grid) &
        &result(ierr)
        use num_vars, only: eq_style, eq_jobs_lims, eq_job_nr
        use grid_vars, only: n_r_eq
        use grid_utilities, only: calc_n_par_X_rich
        use HELENA_vars, only: nchi, chi_H
        use rich_vars, only: n_par_X
        
        character(*), parameter :: rout_name = 'setup_grid_eq'
        
        ! input / output
        type(grid_type), intent(inout) :: grid_eq                               ! equilibrium grid
        integer, intent(in) :: eq_limits(2)                                     ! min. and max. index of eq grid of this process
        logical, intent(in), optional :: only_half_grid                         ! calculate only half grid with even points
        
        ! local variables
        integer :: id                                                           ! counter
        integer :: n_par_X_rich                                                 ! n_par_X for this Richardson level
        
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
                
                ! set local n_par_X
                ierr = calc_n_par_X_rich(n_par_X_rich,only_half_grid)
                CHCKERR('')
                
                ! create grid with eq_jobs_lims
                ierr = grid_eq%init([eq_jobs_lims(2,eq_job_nr)-&
                    &eq_jobs_lims(1,eq_job_nr)+1,1,n_r_eq],eq_limits)           ! only one field line
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
    
    ! Sets up  the field-aligned equilibrium grid,  which serves as a  bridge to
    ! the  solution grid,  as  it contains  the same  normal  coordinate as  the
    ! general  grid, but  the angular  coordinates are  defined by  the solution
    ! grid.
    ! In contrast to setup_grid_eq,  the angular coordinates are also calculated
    ! here.
    integer function setup_grid_eq_B(grid_eq,grid_eq_B,eq,only_half_grid) &
        &result(ierr)
        use num_vars, only: eq_style, eq_jobs_lims, eq_job_nr
        use grid_utilities, only: calc_n_par_X_rich
        
        character(*), parameter :: rout_name = 'setup_grid_eq_B'
        
        ! input / output
        type(grid_type), intent(inout) :: grid_eq                               ! general equilibrium grid
        type(grid_type), intent(inout) :: grid_eq_B                             ! field-aligned equilibrium grid
        type(eq_1_type), intent(in) :: eq                                       ! flux equilibrium variables
        logical, intent(in), optional :: only_half_grid                         ! calculate only half grid with even points
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer :: n_par_X_rich                                                 ! n_par_X for this Richardson level
        
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
                ! set local n_par_X
                ierr = calc_n_par_X_rich(n_par_X_rich,only_half_grid)
                CHCKERR('')
                
                ! create grid with eq_jobs_lims
                ierr = grid_eq_B%init([eq_jobs_lims(2,eq_job_nr)-&
                    &eq_jobs_lims(1,eq_job_nr)+1,1,grid_eq%n(3)],&
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
    
    ! Sets up the general perturbation grid, in which the perturbation variables
    ! are calculated. This  grid has the same angular extent  as the equilibrium
    ! grid but with a higher number  of normal points, indicated by the variable
    ! 'r_F_X'.
    ! Note that no ghost  range is needed as the full normal  range is used: The
    ! division is in the mode numbers
    integer function setup_grid_X(grid_eq,grid_X,r_F_X,X_limits) result(ierr)
        use num_vars, only: norm_disc_prec_X
        use grid_vars, only: disc_type
        use grid_utilities, only: coord_F2E, setup_interp_data, apply_disc
        
        character(*), parameter :: rout_name = 'setup_grid_X'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid
        type(grid_type), intent(inout) :: grid_X                                ! perturbation grid
        real(dp), intent(in) :: r_F_X(:)                                        ! points of perturbation grid
        integer, intent(in) :: X_limits(2)                                      ! min. and max. index of perturbation grid of this process
        
        ! local variables
        type(disc_type) :: norm_interp_data                                     ! data for normal interpolation
        
        ! initialize ierr
        ierr = 0
        
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
        
        ! setup normal interpolation data
        ierr = setup_interp_data(grid_eq%loc_r_F,grid_X%loc_r_F,&
            &norm_interp_data,norm_disc_prec_X)
        CHCKERR('')
        
        ! interpolate
        ierr = apply_disc(grid_eq%theta_E,norm_interp_data,grid_X%theta_E,3)
        CHCKERR('')
        ierr = apply_disc(grid_eq%zeta_E,norm_interp_data,grid_X%zeta_E,3)
        CHCKERR('')
        ierr = apply_disc(grid_eq%theta_F,norm_interp_data,grid_X%theta_F,3)
        CHCKERR('')
        ierr = apply_disc(grid_eq%zeta_F,norm_interp_data,grid_X%zeta_F,3)
        CHCKERR('')
        
        ! clean up
        call norm_interp_data%dealloc()
    end function setup_grid_X
    
    ! Sets  up the general  solution grid, in  which the solution  variables are
    ! calculated.
    integer function setup_grid_sol(grid_eq,grid_sol,r_F_sol,sol_limits) &
        &result(ierr)
        use grid_utilities, only: coord_F2E
        
        character(*), parameter :: rout_name = 'setup_grid_sol'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid
        type(grid_type), intent(inout) :: grid_sol                              ! solution grid
        real(dp), intent(in) :: r_F_sol(:)                                      ! points of solution grid
        integer, intent(in) :: sol_limits(2)                                    ! min. and max. index of sol grid of this process
        
        ! initialize ierr
        ierr = 0
        
        ! create grid
        ierr = grid_sol%init([0,0,size(r_F_sol)],sol_limits)
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
    end function setup_grid_sol
    
    ! Calculate grid that follows magnetic field lines.
    ! Note: The end-points are included for the grids in the parallel direction.
    ! This is to  facilitate working with the trapezoidal rule  or Simpson's 3/8
    ! rule for integration. This is NOT valid in general!
    ! Note: by  setting the flag "only_half_grid",  only the even points  of the
    ! parallel grid are calculated, which is useful for higher Richardson levels
    ! with VMEC so that only new angular  points are calculated and the old ones
    ! reused.
    integer function calc_ang_grid_eq_B(grid_eq,eq,only_half_grid) result(ierr)
        use num_vars, only: use_pol_flux_F, use_pol_flux_E, &
            &eq_style, tol_zero, eq_job_nr, eq_jobs_lims
        use grid_vars, only: min_par_X, max_par_X
        use sol_vars, only: alpha
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
        type(grid_type), intent(inout) :: grid_eq                               ! equilibrium grid of which to calculate angular part
        type(eq_1_type), intent(in), target :: eq                               ! flux equilibrium variables
        logical, intent(in), optional :: only_half_grid                         ! calculate only half grid with even points
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        real(dp), allocatable :: r_E_loc(:)                                     ! flux in Equilibrium coords.
        real(dp), pointer :: flux_F(:) => null()                                ! flux that the F uses as normal coord.
        real(dp), pointer :: flux_E(:) => null()                                ! flux that the E uses as normal coord.
        real(dp) :: r_F_factor, r_E_factor                                      ! mult. factors for r_F and r_E
        integer :: pmone                                                        ! plus or minus one
        integer :: id, kd                                                       ! counters
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
            &range '//trim(r2strt(min_par_X*pi))//'..'//&
            &trim(r2strt(max_par_X*pi)))
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
        allocate(theta_F_loc(n_par_X,1,grid_eq%loc_n_r))
        allocate(zeta_F_loc(n_par_X,1,grid_eq%loc_n_r))
        if (use_pol_flux_F) then                                                ! parallel angle theta
            ierr = calc_eqd_grid(theta_F_loc,min_par_X*pi,&
                &max_par_X*pi,1,excl_last=.false.)                              ! first index corresponds to parallel angle
            CHCKERR('')
            do kd = 1,grid_eq%loc_n_r
                zeta_F_loc(:,:,kd) = pmone*eq%q_saf_E(kd,0)*&
                    &theta_F_loc(:,:,kd)
            end do
            zeta_F_loc = zeta_F_loc + alpha
        else                                                                    ! parallel angle zeta
            ierr = calc_eqd_grid(zeta_F_loc,min_par_X*pi,&
                &max_par_X*pi,1,excl_last=.false.)                              ! first index corresponds to parallel angle
            CHCKERR('')
            do kd = 1,grid_eq%loc_n_r
                theta_F_loc(:,:,kd) = pmone*eq%rot_t_E(kd,0)*&
                    &zeta_F_loc(:,:,kd)
            end do
            theta_F_loc = theta_F_loc - alpha
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
        ierr = coord_F2E(grid_eq,grid_eq%loc_r_F,grid_eq%theta_F,&
            &grid_eq%zeta_F,r_E_loc,grid_eq%theta_E,grid_eq%zeta_E,&
            &r_F_array=flux_F/r_F_factor,r_E_array=flux_E/r_E_factor)
        CHCKERR('')
        
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
    
    ! plots the grid in real space
    ! The  equilibrium grid  should contain  the fieldline-oriented  angles with
    ! ang_1 the parallel angle and ang_2 the field line label.
    ! Note:  This  routine  does  not  use  n_theta_plot  and  n_zeta_plot  from
    ! num_vars, but instead  temporarily overwrites them with its  own, since it
    ! has to be 3D also in the axisymmetric case.
    ! [MPI] Collective call
    integer function magn_grid_plot(grid) result(ierr)
        use num_vars, only: rank, no_plots, n_theta_plot, n_zeta_plot, &
            &eq_style, min_theta_plot, max_theta_plot, min_zeta_plot, &
            &max_zeta_plot
        use grid_utilities, only: trim_grid, extend_grid_E, calc_XYZ_grid
        use VMEC, only: calc_trigon_factors
        
        character(*), parameter :: rout_name = 'magn_grid_plot'
        
        ! input / output
        type(grid_type), intent(in) :: grid                                     ! fieldline-oriented equilibrium grid
        
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
        ierr = extend_grid_E(grid,grid_ext)
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
        call get_full_XYZ(X_1,Y_1,Z_1,X_1_tot,Y_1_tot,Z_1_tot,'flux surfaces')
        call get_full_XYZ(X_2,Y_2,Z_2,X_2_tot,Y_2_tot,Z_2_tot,'field lines')
        
        ierr = magn_grid_plot_HDF5(X_1_tot,X_2_tot,Y_1_tot,Y_2_tot,&
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
        integer function magn_grid_plot_HDF5(X_1,X_2,Y_1,Y_2,Z_1,Z_2,&
            &anim_name) result(ierr)
            use HDF5_ops, only: open_HDF5_file, add_HDF5_item, &
                &print_HDF5_top, print_HDF5_geom, print_HDF5_3D_data_item, &
                &print_HDF5_grid, close_HDF5_file
            use HDF5_vars, only: dealloc_XML_str, &
                & XML_str_type, HDF5_file_type
            use rich_vars, only: rich_lvl
            
            character(*), parameter :: rout_name = 'magn_grid_plot_HDF5'
            
            ! input / output
            real(dp), intent(in) :: X_1(:,:,:), Y_1(:,:,:), Z_1(:,:,:)          ! X, Y and Z of surface in Axisymmetric coordinates
            real(dp), intent(in) :: X_2(:,:,:), Y_2(:,:,:), Z_2(:,:,:)          ! X, Y and Z of magnetic field lines in Axisymmetric coordinates
            character(len=*), intent(in) :: anim_name                           ! name of animation
            
            ! local variables
            character(len=max_str_ln) :: plot_title(2)                          ! plot title for flux surface and for field lines
            character(len=max_str_ln) :: file_name                              ! name of file
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
                
                ! set up plot titles and file name
                plot_title(1) = 'Magnetic Flux Surfaces'
                plot_title(2) = 'Magnetic Field Lines'
                file_name = 'field_lines_in_flux_surfaces_R_'//&
                    &trim(i2str(rich_lvl))
                
                ! open HDF5 file
                ierr = open_HDF5_file(file_info,trim(file_name),&
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
    
    ! Print grid variables to an output file.
    ! Note: "grid_" is added in front the data_name.
    ! If "rich_lvl" is  provided, "_R_rich_lvl" is appended to the  data name if
    ! it is > 0 (only for eq_2), and similarly for "eq_job" through "_E_eq_job".
    integer function print_output_grid(grid,grid_name,data_name,rich_lvl,&
        &eq_job) result(ierr)
        use num_vars, only: PB3D_name, rank
        use HDF5_ops, only: print_HDF5_arrs
        use HDF5_vars, only: var_1D_type, &
            &max_dim_var_1D
        use grid_utilities, only: trim_grid
        
        character(*), parameter :: rout_name = 'print_output_grid'
        
        ! input / output
        type(grid_type), intent(in) :: grid                                     ! grid variables
        character(len=*), intent(in) :: grid_name                               ! name to display
        character(len=*), intent(in) :: data_name                               ! name under which to store
        integer, intent(in), optional :: rich_lvl                               ! Richardson level to reconstruct
        integer, intent(in), optional :: eq_job                                 ! equilibrium job to print
        
        ! local variables
        type(grid_type) :: grid_trim                                            ! trimmed grid
        type(var_1D_type), allocatable, target :: grid_1D(:)                    ! 1D equivalent of grid variables
        type(var_1D_type), pointer :: grid_1D_loc => null()                     ! local element in grid_1D
        integer :: id                                                           ! counter
        
        ! initialize ierr
        ierr = 0
        
        ! user output
        call writo('Write '//trim(grid_name)//' grid variables to output &
            &file')
        call lvl_ud(1)
        
        ! trim grid
        ierr = trim_grid(grid,grid_trim)
        CHCKERR('')
        
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
        if (rank.eq.0) then
            grid_1D_loc%loc_i_min = [1]
            grid_1D_loc%loc_i_max = [3]
            allocate(grid_1D_loc%p(3))
            grid_1D_loc%p = grid_trim%n
        else
            grid_1D_loc%loc_i_min = [1]
            grid_1D_loc%loc_i_max = [0]
            allocate(grid_1D_loc%p(0))
        end if
        
        ! r_F
        grid_1D_loc => grid_1D(id); id = id+1
        grid_1D_loc%var_name = 'r_F'
        allocate(grid_1D_loc%tot_i_min(1),grid_1D_loc%tot_i_max(1))
        allocate(grid_1D_loc%loc_i_min(1),grid_1D_loc%loc_i_max(1))
        grid_1D_loc%tot_i_min = [1]
        grid_1D_loc%tot_i_max = [grid_trim%n(3)]
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
        grid_1D_loc%tot_i_max = [grid_trim%n(3)]
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
            grid_1D_loc%tot_i_max = grid_trim%n
            grid_1D_loc%loc_i_min = [1,1,grid_trim%i_min]
            grid_1D_loc%loc_i_max = [grid_trim%n(1:2),grid_trim%i_max]
            allocate(grid_1D_loc%p(size(grid_trim%theta_F)))
            grid_1D_loc%p = reshape(grid_trim%theta_F,[size(grid_trim%theta_F)])
            
            ! theta_E
            grid_1D_loc => grid_1D(id); id = id+1
            grid_1D_loc%var_name = 'theta_E'
            allocate(grid_1D_loc%tot_i_min(3),grid_1D_loc%tot_i_max(3))
            allocate(grid_1D_loc%loc_i_min(3),grid_1D_loc%loc_i_max(3))
            grid_1D_loc%tot_i_min = [1,1,1]
            grid_1D_loc%tot_i_max = grid_trim%n
            grid_1D_loc%loc_i_min = [1,1,grid_trim%i_min]
            grid_1D_loc%loc_i_max = [grid_trim%n(1:2),grid_trim%i_max]
            allocate(grid_1D_loc%p(size(grid_trim%theta_E)))
            grid_1D_loc%p = reshape(grid_trim%theta_E,[size(grid_trim%theta_E)])
            
            ! zeta_F
            grid_1D_loc => grid_1D(id); id = id+1
            grid_1D_loc%var_name = 'zeta_F'
            allocate(grid_1D_loc%tot_i_min(3),grid_1D_loc%tot_i_max(3))
            allocate(grid_1D_loc%loc_i_min(3),grid_1D_loc%loc_i_max(3))
            grid_1D_loc%tot_i_min = [1,1,1]
            grid_1D_loc%tot_i_max = grid_trim%n
            grid_1D_loc%loc_i_min = [1,1,grid_trim%i_min]
            grid_1D_loc%loc_i_max = [grid_trim%n(1:2),grid_trim%i_max]
            allocate(grid_1D_loc%p(size(grid_trim%zeta_F)))
            grid_1D_loc%p = reshape(grid_trim%zeta_F,[size(grid_trim%zeta_F)])
            
            ! zeta_E
            grid_1D_loc => grid_1D(id); id = id+1
            grid_1D_loc%var_name = 'zeta_E'
            allocate(grid_1D_loc%tot_i_min(3),grid_1D_loc%tot_i_max(3))
            allocate(grid_1D_loc%loc_i_min(3),grid_1D_loc%loc_i_max(3))
            grid_1D_loc%tot_i_min = [1,1,1]
            grid_1D_loc%tot_i_max = grid_trim%n
            grid_1D_loc%loc_i_min = [1,1,grid_trim%i_min]
            grid_1D_loc%loc_i_max = [grid_trim%n(1:2),grid_trim%i_max]
            allocate(grid_1D_loc%p(size(grid_trim%zeta_E)))
            grid_1D_loc%p = reshape(grid_trim%zeta_E,[size(grid_trim%zeta_E)])
        end if
        
        ! write
        ierr = print_HDF5_arrs(grid_1D(1:id-1),PB3D_name,&
            &'grid_'//trim(data_name),rich_lvl=rich_lvl,eq_job=eq_job)
        CHCKERR('')
        
        ! clean up
        call grid_trim%dealloc()
        
        ! clean up
        nullify(grid_1D_loc)
        
        ! user output
        call lvl_ud(-1)
    end function print_output_grid
end module grid_ops
