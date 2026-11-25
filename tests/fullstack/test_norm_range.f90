!------------------------------------------------------------------------------!
!> Full-stack tests of the end points of the solution grid (issue #88).
!!
!! calc_norm_range() selects the equilibrium points of a run with con2dis(),
!! which takes a limit within a relative tol_con2dis of an equilibrium point
!! for that point. When min_r_sol or max_r_sol is the decimal value of an
!! equilibrium point, the first or last point of the solution grid can lie
!! outside the equilibrium grid by rounding. setup_grid_sol() (X_grid_style 1
!! and 3) and setup_grid_X() (X_grid_style 2) put such a point on the end
!! point of the equilibrium grid. The tests:
!!
!!  - min_r_sol on each point of the cbm18a fixture with 0.60 <= s <= 0.80,
!!    and max_r_sol on each point with 0.80 <= s <= 0.90, s being the square
!!    root of the normalized poloidal flux, with the flux of the file and with
!!    the flux negated (a decreasing normal coordinate), for X_grid_style 1, 2
!!    and 3: the grids are set up, the end points of the solution grid that
!!    lay outside the equilibrium grid are on its end points, and all other
!!    points are unchanged. At least one case has to have an end point
!!    outside, else the test does not test the rounding;
!!  - a synthetic equilibrium grid, increasing and decreasing, and a solution
!!    grid with its first or last point outside it by one unit in the last
!!    place, by half of tol_con2dis and by twice tol_con2dis of the end point,
!!    and by 1e-6 of it: the first two are put on the end point, the last two
!!    make the set-up return an error.
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
            &new_unittest("tolerance", test_tolerance) &
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

    !> Sets min_r_sol (\c on_min) or max_r_sol to \f$(k/800)^2\f$, the
    !! normalized poloidal flux of point k+1 of the fixture, for k = \c k_lo ..
    !! \c k_hi, and the other limit to \c other_lim, and sets up the
    !! solution grid as PB3D does.
    subroutine scan_points(error,on_min,k_lo,k_hi,other_lim)
        use num_vars, only: eq_i, eq_name, eq_style, use_pol_flux_E, &
            &use_pol_flux_F, norm_disc_prec_eq, norm_disc_prec_X, &
            &norm_disc_prec_sol, sol_n_procs, n_procs, X_grid_style, no_output
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
        integer :: n_out(2)                                                     ! nr. of cases with an end point outside, per sign
        integer :: in_limits(2), sol_limits(2)                                  ! limits of equilibrium and solution grid
        character(len=1024) :: fix_dir                                          ! fixture directory
        character(len=200) :: fail_msg                                          ! first failure
        logical :: no_output_saved                                              ! no_output on entry
        real(dp) :: lim                                                         ! limit on point k+1
        real(dp) :: r_eq_lim(2)                                                 ! lower and upper end of the equilibrium grid
        real(dp), allocatable :: r_F_sol(:)                                     ! solution grid of calc_norm_range
        real(dp), allocatable :: r_F_X(:)                                       ! perturbation grid for X_grid_style 2
        real(dp), allocatable :: r_F_exp(:)                                     ! expected solution grid after the set-up
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
        sol_n_procs = n_procs

        no_output_saved = no_output
        no_output = .true.
        fail_msg = ''
        n_out = 0
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

                ! equilibrium grid of the run
                ierr = calc_norm_range('PB3D_in',in_limits=in_limits)
                if (ierr.ne.0) then
                    fail_msg = 'calc_norm_range PB3D_in'
                    exit signs
                end if
                n_r_eq = in_limits(2)-in_limits(1)+1
                ierr = grid_eq%init([1,1,n_r_eq])
                if (ierr.ne.0) then
                    fail_msg = 'init of grid_eq'
                    exit signs
                end if
                grid_eq%r_F = flux_p_H(in_limits(1):in_limits(2),0)/(2*pi)     ! as calc_flux_q_HEL
                grid_eq%r_E = grid_eq%r_F
                grid_eq%theta_E = 0._dp
                grid_eq%zeta_E = 0._dp
                grid_eq%theta_F = 0._dp
                grid_eq%zeta_F = 0._dp
                r_eq_lim = [minval(grid_eq%r_F),maxval(grid_eq%r_F)]

                ! solution grid of the run
                allocate(r_F_sol(n_r_sol))
                ierr = calc_norm_range('PB3D_sol',sol_limits=sol_limits,&
                    &r_F_sol=r_F_sol)
                if (ierr.ne.0) then
                    fail_msg = 'calc_norm_range PB3D_sol'
                    exit signs
                end if
                if (minval(r_F_sol).lt.r_eq_lim(1) .or. &
                    &maxval(r_F_sol).gt.r_eq_lim(2)) &
                    &n_out(i_sign) = n_out(i_sign) + 1
                r_F_exp = min(max(r_F_sol,r_eq_lim(1)),r_eq_lim(2))

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
                    if (any(grid_sol%r_F.ne.r_F_exp) .or. &
                        &any(grid_sol%loc_r_F.ne.&
                        &r_F_exp(sol_limits(1):sol_limits(2)))) then
                        fail_msg = 'solution grid not as expected'
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
                &', X_grid_style = '//trim(i2str(X_grid_style))//&
                &', flux sign '//trim(i2str(3-2*i_sign))
        end if
        call check(error, len_trim(fail_msg), 0, trim(fail_msg))
        if (allocated(error)) return
        write(error_unit,'(A,I0,A,I0,A)') '   [norm_range] ',n_out(1),' of ',&
            &k_hi-k_lo+1,' cases with an end point outside before the set-up'
        call check(error, n_out(1).gt.0, &
            &'no case with an end point outside the equilibrium grid')
        if (allocated(error)) return
        call check(error, n_out(2), n_out(1), &
            &'negated flux: other number of cases outside')
        if (allocated(error)) return
    end subroutine scan_points

    !> A synthetic equilibrium grid of 9 points, increasing and decreasing, and
    !! a solution grid of 7 points with its first or last point outside it.
    subroutine test_tolerance(error)
        use num_vars, only: n_procs, norm_disc_prec_eq, norm_disc_prec_X, &
            &X_grid_style, no_output
        use num_utilities, only: tol_con2dis
        use grid_ops, only: setup_grid_X, setup_grid_sol

        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: n_eq = 9                                          ! nr. of equilibrium points
        integer, parameter :: n_sol = 7                                         ! nr. of solution points
        integer :: ierr                                                         ! error status
        integer :: i_sign, i_end, i_case, i_style, id                           ! counters
        integer :: ie, is                                                       ! index of the end point in eq. and sol. grid
        character(len=200) :: fail_msg                                          ! first failure
        logical :: no_output_saved                                              ! no_output on entry
        logical :: accept                                                       ! the case has to be accepted
        real(dp) :: dir                                                         ! outward direction at the end point
        real(dp) :: r_eq(n_eq)                                                  ! equilibrium grid
        real(dp), allocatable :: r_sol(:)                                       ! solution grid
        real(dp), allocatable :: r_exp(:)                                       ! expected solution grid after the set-up
        type(grid_type) :: grid_eq, grid_X, grid_sol                            ! grids

        norm_disc_prec_eq = 3
        norm_disc_prec_X = 3

        no_output_saved = no_output
        no_output = .true.
        fail_msg = ''
        cases: do i_sign = 1,2
            r_eq = (3-2*i_sign)*[(0.3_dp+0.05_dp*(id-1),id=1,n_eq)]
            do i_end = 1,2
                if (i_end.eq.1) then
                    ie = 1
                    is = 1
                    dir = sign(1._dp,r_eq(1)-r_eq(n_eq))
                else
                    ie = n_eq
                    is = n_sol
                    dir = sign(1._dp,r_eq(n_eq)-r_eq(1))
                end if
                do i_case = 1,4
                    r_sol = [(r_eq(1)+(id-1)*(r_eq(n_eq)-r_eq(1))/(n_sol-1),&
                        &id=1,n_sol)]
                    r_sol([1,n_sol]) = r_eq([1,n_eq])
                    r_exp = r_sol
                    r_exp(is) = r_eq(ie)
                    select case (i_case)
                        case (1)                                                ! one unit in the last place
                            r_sol(is) = nearest(r_eq(ie),dir)
                        case (2)                                                ! half the tolerance
                            r_sol(is) = r_eq(ie) + &
                                &dir*0.5_dp*tol_con2dis*abs(r_eq(ie))
                        case (3)                                                ! twice the tolerance
                            r_sol(is) = r_eq(ie) + &
                                &dir*2._dp*tol_con2dis*abs(r_eq(ie))
                        case (4)                                                ! 1e-6
                            r_sol(is) = r_eq(ie) + dir*1.E-6_dp*abs(r_eq(ie))
                    end select
                    accept = i_case.le.2

                    do i_style = 1,2
                        X_grid_style = i_style
                        ierr = grid_eq%init([1,1,n_eq])
                        if (ierr.ne.0) then
                            fail_msg = 'init of grid_eq'
                            exit cases
                        end if
                        grid_eq%r_F = r_eq
                        grid_eq%r_E = 2*r_eq
                        grid_eq%theta_E = 0._dp
                        grid_eq%zeta_E = 0._dp
                        grid_eq%theta_F = 0._dp
                        grid_eq%zeta_F = 0._dp
                        if (X_grid_style.eq.2) then
                            ierr = setup_grid_X(grid_eq,grid_X,r_sol,&
                                &[1,n_sol])
                            if (ierr.eq.0) ierr = setup_grid_sol(grid_eq,&
                                &grid_X,grid_sol,r_sol,[1,n_sol])
                            call grid_X%dealloc()
                        else
                            ierr = setup_grid_sol(grid_eq,grid_X,grid_sol,&
                                &r_sol,[1,n_sol])
                        end if
                        if (accept) then
                            if (ierr.ne.0) then
                                fail_msg = 'point within the tolerance refused'
                                exit cases
                            end if
                            if (any(grid_sol%r_F.ne.r_exp) .or. &
                                &abs(grid_sol%r_E(is)-2*r_eq(ie)).gt.&
                                &1.E-14_dp*abs(r_eq(ie))) then
                                fail_msg = 'point not put on the end point'
                                exit cases
                            end if
                        else if (ierr.eq.0) then
                            fail_msg = 'point outside the tolerance accepted'
                            exit cases
                        end if
                        if (associated(grid_sol%r_F)) call grid_sol%dealloc()
                        call grid_eq%dealloc()
                    end do
                end do
            end do
        end do cases
        no_output = no_output_saved

        if (len_trim(fail_msg).gt.0) then
            fail_msg = trim(fail_msg)//': flux sign '//&
                &trim(i2str(3-2*i_sign))//', end '//trim(i2str(i_end))//&
                &', case '//trim(i2str(i_case))//', X_grid_style '//&
                &trim(i2str(X_grid_style))
        end if
        call check(error, len_trim(fail_msg), 0, trim(fail_msg))
    end subroutine test_tolerance
end module test_norm_range
