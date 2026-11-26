!------------------------------------------------------------------------------!
!> Full-stack tests of the equilibrium range of a run (issue #88).
!!
!! calc_norm_range() chooses the equilibrium points of a run so that they
!! contain the solution grid: from the last point at or below its first point
!! to the first point at or above its last point, compared in the numbers
!! that setup_grid_sol() and setup_grid_X() later compare. The tests:
!!
!!  - min_r_sol on each point of the cbm18a fixture with 0.60 <= s <= 0.80,
!!    and max_r_sol on each point with 0.80 <= s <= 0.90, s being the square
!!    root of the normalized poloidal flux, with the flux of the file and with
!!    the flux negated (a decreasing normal coordinate), for X_grid_style 1, 2
!!    and 3: the equilibrium grid contains the solution grid and has no point
!!    to spare at either end, the grids are set up, and the solution grid is
!!    the one of calc_norm_range(). At least one case has to take the point
!!    beyond the limit's point, else the test does not test the rounding;
!!  - a synthetic flux, increasing and decreasing, with min_r_sol or
!!    max_r_sol on each inner point, one unit in the last place on either
!!    side of it, and 1e-9 of its value on either side: the same conditions,
!!    for X_grid_style 1.
!!
!! The suite sets the global input variables that these routines read. The
!! fixture directory is passed through the environment variable
!! PB3D_FIXTURE_DIR (set by CTest; the suite is only registered when the
!! fixture is available).
!------------------------------------------------------------------------------!
module test_norm_range
    use, intrinsic :: iso_fortran_env, only: error_unit
    use testdrive, only: new_unittest, unittest_type, error_type, check
    use num_vars, only: dp, pi
    use str_utilities, only: i2str
    use grid_vars, only: grid_type

    implicit none
    private
    public collect_norm_range

    integer, parameter :: n_s = 800                                             ! the fixture has its points at s = 0, 1/n_s, .., 1

contains
    !> Collect the tests.
    subroutine collect_norm_range(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
            &new_unittest("min_r_sol_on_points", test_min_on_points), &
            &new_unittest("max_r_sol_on_points", test_max_on_points), &
            &new_unittest("limits_near_points", test_near_points) &
            &]
    end subroutine collect_norm_range

    !> min_r_sol on the points with s = 0.60 .. 0.80, max_r_sol = 0.75.
    subroutine test_min_on_points(error)
        type(error_type), allocatable, intent(out) :: error

        call scan_points(error,.true.,480,640,0.75_dp)
    end subroutine test_min_on_points

    !> max_r_sol on the points with s = 0.80 .. 0.90, min_r_sol = 0.45.
    subroutine test_max_on_points(error)
        type(error_type), allocatable, intent(out) :: error

        call scan_points(error,.false.,640,720,0.45_dp)
    end subroutine test_max_on_points

    !> Checks the equilibrium range \c in_limits of the flux \c flux against
    !! the solution grid \c r_F_sol: the equilibrium grid flux/2pi contains
    !! the solution grid, and removing its first or last point would leave a
    !! point of the solution grid outside it. Returns an empty message if
    !! both hold.
    function range_msg(flux,in_limits,r_F_sol) result(msg)
        real(dp), intent(in) :: flux(:)                                         !< flux of all equilibrium points
        integer, intent(in) :: in_limits(2)                                     !< equilibrium range
        real(dp), intent(in) :: r_F_sol(:)                                      !< solution grid
        character(len=100) :: msg

        real(dp) :: pm                                                          ! sign of the flux
        real(dp), allocatable :: r_eq(:)                                        ! pm flux/2pi of the equilibrium range
        real(dp) :: r_sol(2)                                                    ! pm times first and last solution point

        msg = ''
        pm = sign(1._dp,flux(size(flux)))
        r_eq = pm*(flux(in_limits(1):in_limits(2))/(2*pi))
        r_sol = pm*r_F_sol([1,size(r_F_sol)])
        if (r_eq(1).gt.r_sol(1) .or. r_eq(size(r_eq)).lt.r_sol(2)) then
            msg = 'solution grid outside the equilibrium grid'
        else if (size(r_eq).gt.1) then
            if (r_eq(2).le.r_sol(1) .or. r_eq(size(r_eq)-1).ge.r_sol(2)) &
                &msg = 'equilibrium grid larger than needed'
        end if
    end function range_msg

    !> Sets min_r_sol (\c on_min) or max_r_sol to \f$(k/800)^2\f$, the
    !! normalized poloidal flux of point k+1 of the fixture, for k = \c k_lo ..
    !! \c k_hi, and the other limit to \c other_lim, and sets up the
    !! solution grid as PB3D does.
    subroutine scan_points(error,on_min,k_lo,k_hi,other_lim)
        use num_vars, only: eq_i, eq_name, eq_style, use_pol_flux_E, &
            &use_pol_flux_F, norm_disc_prec_eq, norm_disc_prec_X, &
            &norm_disc_prec_sol, X_grid_style, no_output
        use grid_vars, only: n_r_in, n_r_eq, n_r_sol
        use X_vars, only: min_r_sol, max_r_sol
        use HELENA_vars, only: flux_p_H, dealloc_HEL
        use HELENA_ops, only: read_HEL
        use grid_ops, only: calc_norm_range, setup_grid_X, setup_grid_sol

        type(error_type), allocatable, intent(out) :: error
        logical, intent(in) :: on_min                                           !< set min_r_sol, else max_r_sol
        integer, intent(in) :: k_lo, k_hi                                       !< range of points
        real(dp), intent(in) :: other_lim                                       !< value of the other limit

        integer :: ierr, istat                                                  ! error status
        integer :: i_sign, k, i_style                                           ! counters
        integer :: n_beyond(2)                                                  ! cases beyond the limit's point, per sign
        integer :: in_limits(2), sol_limits(2)                                  ! limits of equilibrium and solution grid
        character(len=1024) :: fix_dir                                          ! fixture directory
        character(len=200) :: fail_msg                                          ! first failure
        logical :: no_output_saved                                              ! no_output on entry
        real(dp) :: lim                                                         ! limit on point k+1
        real(dp), allocatable :: r_F_sol(:)                                     ! solution grid of calc_norm_range
        real(dp), allocatable :: r_F_X(:)                                       ! perturbation grid for X_grid_style 2
        type(grid_type) :: grid_eq, grid_X, grid_sol                            ! grids

        ! read the fixture
        call get_environment_variable('PB3D_FIXTURE_DIR',fix_dir,status=istat)
        call check(error, istat, 0, 'PB3D_FIXTURE_DIR not set')
        if (allocated(error)) return
        eq_name = trim(fix_dir)//'/cbm18a_extended_vac.12'
        open(unit=eq_i,file=trim(eq_name),status='old',action='read',&
            &iostat=istat)
        call check(error, istat, 0, 'could not open '//trim(eq_name))
        if (allocated(error)) return
        if (allocated(flux_p_H)) call dealloc_HEL
        ierr = read_HEL(n_r_in,use_pol_flux_E)
        close(eq_i)
        call check(error, ierr, 0, 'read_HEL failed')
        if (allocated(error)) return
        call check(error, n_r_in, n_s+1, 'number of points of the fixture')
        if (allocated(error)) return

        ! input as in the regression deck input_cbm18a_fixed_bnd
        eq_style = 2
        use_pol_flux_F = .true.
        norm_disc_prec_eq = 3
        norm_disc_prec_X = 3
        norm_disc_prec_sol = 3
        n_r_sol = 200

        no_output_saved = no_output
        no_output = .true.
        fail_msg = ''
        n_beyond = 0
        signs: do i_sign = 1,2
            if (i_sign.eq.2) flux_p_H = -flux_p_H
            do k = k_lo,k_hi
                lim = real(k,dp)**2/real(n_s,dp)**2
                if (on_min) then
                    min_r_sol = lim
                    max_r_sol = other_lim
                else
                    min_r_sol = other_lim
                    max_r_sol = lim
                end if

                ! equilibrium range and solution grid of the run
                i_style = 0
                ierr = calc_norm_range('PB3D_in',in_limits=in_limits)
                if (ierr.ne.0) then
                    fail_msg = 'calc_norm_range PB3D_in'
                    exit signs
                end if
                allocate(r_F_sol(n_r_sol))
                ierr = calc_norm_range('PB3D_sol',sol_limits=sol_limits,&
                    &r_F_sol=r_F_sol)
                if (ierr.ne.0) then
                    fail_msg = 'calc_norm_range PB3D_sol'
                    exit signs
                end if
                fail_msg = range_msg(flux_p_H(:,0),in_limits,r_F_sol)
                if (len_trim(fail_msg).gt.0) exit signs
                if (on_min .and. in_limits(1).eq.k .or. &
                    &.not.on_min .and. in_limits(2).eq.k+2) &
                    &n_beyond(i_sign) = n_beyond(i_sign) + 1

                ! equilibrium grid as calc_flux_q_HEL sets it up
                n_r_eq = in_limits(2)-in_limits(1)+1
                ierr = grid_eq%init([1,1,n_r_eq],divided=.false.)
                if (ierr.ne.0) then
                    fail_msg = 'init of grid_eq'
                    exit signs
                end if
                grid_eq%r_F = flux_p_H(in_limits(1):in_limits(2),0)/(2*pi)
                grid_eq%r_E = grid_eq%r_F
                grid_eq%theta_E = 0._dp
                grid_eq%zeta_E = 0._dp
                grid_eq%theta_F = 0._dp
                grid_eq%zeta_F = 0._dp

                do i_style = 1,3
                    X_grid_style = i_style
                    if (X_grid_style.eq.2) then
                        r_F_X = r_F_sol
                        ierr = setup_grid_X(grid_eq,grid_X,r_F_X,sol_limits)
                        if (ierr.ne.0) then
                            fail_msg = 'setup_grid_X'
                            exit signs
                        end if
                    end if
                    ierr = setup_grid_sol(grid_eq,grid_X,grid_sol,r_F_sol,&
                        &sol_limits)
                    if (ierr.ne.0) then
                        fail_msg = 'setup_grid_sol'
                        exit signs
                    end if
                    if (any(grid_sol%r_F.ne.r_F_sol)) then
                        fail_msg = 'solution grid changed'
                        exit signs
                    end if
                    call grid_sol%dealloc()
                    if (X_grid_style.eq.2) call grid_X%dealloc()
                end do
                call grid_eq%dealloc()
                deallocate(r_F_sol)
            end do
        end do signs
        no_output = no_output_saved
        call dealloc_HEL

        if (len_trim(fail_msg).gt.0) then
            fail_msg = trim(fail_msg)//', k = '//trim(i2str(k))//&
                &', flux sign '//trim(i2str(3-2*i_sign))
            if (i_style.gt.0) fail_msg = trim(fail_msg)//&
                &', X_grid_style = '//trim(i2str(i_style))
        end if
        call check(error, len_trim(fail_msg), 0, trim(fail_msg))
        if (allocated(error)) return
        write(error_unit,'(A,I0,A,I0,A)') '   [norm_range] ',n_beyond(1),&
            &' of ',k_hi-k_lo+1,' cases take the point beyond the limit''s point'
        call check(error, n_beyond(1).gt.0, &
            &'no case takes the point beyond the limit''s point')
        if (allocated(error)) return
        call check(error, n_beyond(2), n_beyond(1), &
            &'negated flux: other number of cases beyond')
        if (allocated(error)) return
    end subroutine scan_points

    !> A synthetic flux of 21 points, increasing and decreasing, with
    !! min_r_sol or max_r_sol on each inner point, one unit in the last place
    !! on either side of it, and 1e-9 of its value on either side.
    subroutine test_near_points(error)
        use num_vars, only: eq_style, use_pol_flux_E, use_pol_flux_F, &
            &norm_disc_prec_eq, norm_disc_prec_sol, X_grid_style, no_output
        use grid_vars, only: n_r_in, n_r_eq, n_r_sol
        use X_vars, only: min_r_sol, max_r_sol
        use HELENA_vars, only: flux_p_H, flux_t_H
        use grid_ops, only: calc_norm_range, setup_grid_sol

        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: n = 21                                            ! nr. of equilibrium points
        integer :: ierr                                                         ! error status
        integer :: i_sign, i_end, k, i_off, id                                  ! counters
        integer :: in_limits(2), sol_limits(2)                                  ! limits of equilibrium and solution grid
        character(len=200) :: fail_msg                                          ! first failure
        logical :: no_output_saved                                              ! no_output on entry
        real(dp) :: v, lim                                                      ! normalized flux of point k, limit
        real(dp), allocatable :: r_F_sol(:)                                     ! solution grid
        type(grid_type) :: grid_eq, grid_X, grid_sol                            ! grids

        eq_style = 2
        use_pol_flux_E = .true.
        use_pol_flux_F = .true.
        norm_disc_prec_eq = 3
        norm_disc_prec_sol = 3
        X_grid_style = 1
        n_r_in = n
        n_r_sol = 7
        if (allocated(flux_p_H)) deallocate(flux_p_H)
        if (allocated(flux_t_H)) deallocate(flux_t_H)
        allocate(flux_p_H(n,0:0),flux_t_H(n,0:0))

        no_output_saved = no_output
        no_output = .true.
        fail_msg = ''
        cases: do i_sign = 1,2
            flux_p_H(:,0) = (3-2*i_sign)*0.7_dp*[(((id-1._dp)/(n-1))**1.5_dp,&
                &id=1,n)]
            flux_t_H = 1.1_dp*flux_p_H
            do i_end = 1,2
                do k = 2,n-1
                    v = flux_p_H(k,0)/flux_p_H(n,0)
                    do i_off = 1,5
                        select case (i_off)
                            case (1)
                                lim = v
                            case (2)
                                lim = nearest(v,-1._dp)
                            case (3)
                                lim = nearest(v,1._dp)
                            case (4)
                                lim = v*(1-1.E-9_dp)
                            case (5)
                                lim = v*(1+1.E-9_dp)
                        end select
                        if (i_end.eq.1) then
                            min_r_sol = lim
                            max_r_sol = 1._dp
                        else
                            min_r_sol = 0._dp
                            max_r_sol = lim
                        end if

                        ierr = calc_norm_range('PB3D_in',in_limits=in_limits)
                        if (ierr.ne.0) then
                            fail_msg = 'calc_norm_range PB3D_in'
                            exit cases
                        end if
                        allocate(r_F_sol(n_r_sol))
                        ierr = calc_norm_range('PB3D_sol',&
                            &sol_limits=sol_limits,r_F_sol=r_F_sol)
                        if (ierr.ne.0) then
                            fail_msg = 'calc_norm_range PB3D_sol'
                            exit cases
                        end if
                        fail_msg = range_msg(flux_p_H(:,0),in_limits,r_F_sol)
                        if (len_trim(fail_msg).gt.0) exit cases

                        n_r_eq = in_limits(2)-in_limits(1)+1
                        ierr = grid_eq%init([1,1,n_r_eq],divided=.false.)
                        if (ierr.ne.0) then
                            fail_msg = 'init of grid_eq'
                            exit cases
                        end if
                        grid_eq%r_F = flux_p_H(in_limits(1):in_limits(2),0)/&
                            &(2*pi)
                        grid_eq%r_E = grid_eq%r_F
                        grid_eq%theta_E = 0._dp
                        grid_eq%zeta_E = 0._dp
                        grid_eq%theta_F = 0._dp
                        grid_eq%zeta_F = 0._dp
                        ierr = setup_grid_sol(grid_eq,grid_X,grid_sol,r_F_sol,&
                            &sol_limits)
                        if (ierr.ne.0) then
                            fail_msg = 'setup_grid_sol'
                            exit cases
                        end if
                        call grid_sol%dealloc()
                        call grid_eq%dealloc()
                        deallocate(r_F_sol)
                    end do
                end do
            end do
        end do cases
        no_output = no_output_saved
        deallocate(flux_p_H,flux_t_H)

        if (len_trim(fail_msg).gt.0) then
            fail_msg = trim(fail_msg)//': flux sign '//&
                &trim(i2str(3-2*i_sign))//', end '//trim(i2str(i_end))//&
                &', point '//trim(i2str(k))//', offset '//trim(i2str(i_off))
        end if
        call check(error, len_trim(fail_msg), 0, trim(fail_msg))
    end subroutine test_near_points
end module test_norm_range
