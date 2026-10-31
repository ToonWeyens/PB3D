!------------------------------------------------------------------------------!
!> Golden-file test of the HELENA equilibrium parser (HELENA_ops.read_HEL)
!! on the committed cbm18a_extended_vac fixture (the same file the physics
!! regression layer runs the full code on).
!!
!! read_HEL is deterministic given the file, so the parsed and derived
!! quantities (grid sizes, profiles, MISHKA-normalization factors, geometry)
!! are pinned against golden values recorded from the current
!! implementation; internal-consistency relations that must hold exactly by
!! construction (rot_t = 1/q, flux_t' = q flux_p', flux_p' = 2 pi) and an
!! implementation-independent trapezoidal cross-check of the toroidal flux
!! integral guard against parsing- and normalization regressions with
!! independent references.
!!
!! The fixture directory is passed through the environment variable
!! PB3D_FIXTURE_DIR (set by CTest; the test is only registered when the
!! fixture is available).
!------------------------------------------------------------------------------!
module test_read_HEL
    use, intrinsic :: iso_fortran_env, only: error_unit
    use testdrive, only: new_unittest, unittest_type, error_type, check
    use num_vars, only: dp, pi
    use str_utilities, only: r2str, i2str

    implicit none
    private
    public collect_read_HEL

contains
    !> Collect the tests.
    subroutine collect_read_HEL(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
            &new_unittest("cbm18a_golden", test_cbm18a_golden) &
            &]
    end subroutine collect_read_HEL

    !> Parse the fixture and pin everything.
    subroutine test_cbm18a_golden(error)
        use num_vars, only: eq_i, eq_name
        use HELENA_ops, only: read_HEL
        use HELENA_vars, only: pres_H, q_saf_H, rot_t_H, flux_p_H, flux_t_H, &
            &nchi, chi_H, ias, RBphi_H, R_H, Z_H, RMtoG_H, BMtoG_H

        type(error_type), allocatable, intent(out) :: error

        integer :: ierr, kd                                                     ! error status, counter
        integer :: n_r_in                                                       ! number of normal points
        integer :: istat                                                        ! status
        logical :: use_pol_flux_H                                               ! whether HELENA uses poloidal flux
        character(len=1024) :: fix_dir                                          ! fixture directory
        real(dp) :: flux_t_trap                                                 ! trapezoidal toroidal flux
        real(dp) :: ellip                                                       ! ellipticity

        ! open the fixture on the equilibrium unit
        call get_environment_variable('PB3D_FIXTURE_DIR',fix_dir,&
            &status=istat)
        call check(error, istat, 0, 'PB3D_FIXTURE_DIR not set')
        if (allocated(error)) return
        eq_name = trim(fix_dir)//'/cbm18a_extended_vac.12'
        open(unit=eq_i,file=trim(eq_name),status='old',action='read',&
            &iostat=istat)
        call check(error, istat, 0, 'could not open '//trim(eq_name))
        if (allocated(error)) return

        ierr = read_HEL(n_r_in,use_pol_flux_H)
        close(eq_i)
        call check(error, ierr, 0, 'read_HEL failed')
        if (allocated(error)) return

        write(error_unit,'(A,I5,A,I5,A,I2)') '   [read_HEL] n_r_in = ',&
            &n_r_in,', nchi = ',nchi,', ias = ',ias
        write(error_unit,'(A,2ES23.15)') '   [read_HEL] q axis/edge:    ',&
            &q_saf_H(1,0),q_saf_H(n_r_in,0)
        write(error_unit,'(A,2ES23.15)') '   [read_HEL] pres axis/edge: ',&
            &pres_H(1,0),pres_H(n_r_in,0)
        write(error_unit,'(A,2ES23.15)') '   [read_HEL] F axis/edge:    ',&
            &RBphi_H(1,0),RBphi_H(n_r_in,0)
        write(error_unit,'(A,2ES23.15)') '   [read_HEL] flux_p/t edge:  ',&
            &flux_p_H(n_r_in,0),flux_t_H(n_r_in,0)
        write(error_unit,'(A,2ES23.15)') '   [read_HEL] RMtoG, BMtoG:   ',&
            &RMtoG_H,BMtoG_H
        write(error_unit,'(A,3ES23.15)') '   [read_HEL] R lims, max Z:  ',&
            &minval(R_H),maxval(R_H),maxval(Z_H)

        ! grid and symmetry metadata (from the file header)
        call check(error, ias, 0, 'cbm18a is up-down symmetric')
        if (allocated(error)) return
        call check(error, n_r_in, 801, 'n_r_in golden')
        if (allocated(error)) return
        call check(error, nchi, 401, 'nchi golden')
        if (allocated(error)) return
        call check(error, use_pol_flux_H, 'HELENA uses poloidal flux')
        if (allocated(error)) return

        ! symmetric equilibria are given on the half poloidal circumference
        ! (the file stores the angles with 8 significant digits)
        call check(error, abs(chi_H(1)).lt.1.e-12_dp .and. &
            &abs(chi_H(nchi)-pi).lt.1.e-6_dp, &
            &'chi grid spans [0, pi] for ias = 0')
        if (allocated(error)) return

        ! exact-by-construction consistency relations
        call check(error, maxval(abs(flux_p_H(:,1)-2._dp*pi)).lt.1.e-12_dp, &
            &'flux_p_H'' = 2 pi (MISHKA normalization)')
        if (allocated(error)) return
        call check(error, maxval(abs(flux_t_H(:,1)-&
            &q_saf_H(:,0)*2._dp*pi)).lt.1.e-10_dp, &
            &'flux_t'' = q flux_p''')
        if (allocated(error)) return
        call check(error, maxval(abs(rot_t_H(:,0)*q_saf_H(:,0)-1._dp)).lt.&
            &1.e-12_dp, 'rot_t = 1/q')
        if (allocated(error)) return

        ! monotonicity of the fluxes
        do kd = 2,n_r_in
            call check(error, flux_p_H(kd,0).gt.flux_p_H(kd-1,0) .and. &
                &flux_t_H(kd,0).gt.flux_t_H(kd-1,0), &
                &'fluxes increase monotonically')
            if (allocated(error)) return
        end do

        ! implementation-independent cross-check: the toroidal flux is the
        ! integral of q over the poloidal flux (trapezoid, so only accurate
        ! to the grid resolution)
        flux_t_trap = 0._dp
        do kd = 2,n_r_in
            flux_t_trap = flux_t_trap + &
                &0.5_dp*(q_saf_H(kd,0)+q_saf_H(kd-1,0))*&
                &(flux_p_H(kd,0)-flux_p_H(kd-1,0))
        end do
        call check(error, abs(flux_t_H(n_r_in,0)-flux_t_trap).lt.&
            &1.e-3_dp*abs(flux_t_trap), &
            &'flux_t edge vs trapezoidal integral of q dpsi: '//&
            &trim(r2str(flux_t_H(n_r_in,0)))//' vs '//&
            &trim(r2str(flux_t_trap)))
        if (allocated(error)) return

        ! geometry: degenerate first normal point on the magnetic axis
        call check(error, maxval(abs(R_H(:,1)-R_H(1,1))).lt.1.e-12_dp .and. &
            &maxval(abs(Z_H(:,1))).lt.1.e-12_dp, &
            &'first normal point degenerates to the magnetic axis')
        if (allocated(error)) return

        ! golden profile and normalization values (recorded from the
        ! current implementation; Linux gfortran 13, 2026-07-10 - if these
        ! move, the parser or the MISHKA normalization changed)
        call check(error, q_saf_H(1,0), 1.0521365_dp, thr=1.e-10_dp, &
            &rel=.true., message='q on axis golden')
        if (allocated(error)) return
        call check(error, q_saf_H(n_r_in,0), 25.031802_dp, thr=1.e-10_dp, &
            &rel=.true., message='q at edge golden')
        if (allocated(error)) return
        call check(error, pres_H(1,0), 1.2251493e-2_dp, thr=1.e-10_dp, &
            &rel=.true., message='pressure on axis golden')
        if (allocated(error)) return
        call check(error, flux_p_H(n_r_in,0), &
            &3.767600982734008e-1_dp, thr=1.e-10_dp, rel=.true., &
            &message='poloidal flux at edge golden')
        if (allocated(error)) return
        call check(error, RMtoG_H, 8.319203019203020e-1_dp, thr=1.e-10_dp, &
            &rel=.true., message='RMtoG golden')
        if (allocated(error)) return
        call check(error, BMtoG_H, 1.1634848_dp, thr=1.e-10_dp, rel=.true., &
            &message='BMtoG golden')
        if (allocated(error)) return

        ! ellipticity from the parsed flux surfaces (ias = 0: top half)
        ellip = 2._dp*maxval(Z_H)/(maxval(R_H)-minval(R_H))
        call check(error, ellip, 9.9995911e-1_dp, thr=1.e-8_dp, rel=.true., &
            &message='ellipticity golden')
        if (allocated(error)) return
    end subroutine test_cbm18a_golden
end module test_read_HEL
