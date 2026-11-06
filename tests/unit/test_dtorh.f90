!------------------------------------------------------------------------------!
!> Unit tests for dtorh: toroidal harmonics
!! \f$P_{n-1/2}(z)\f$ and \f$Q_{n-1/2}(z)\f$ (order m = 0).
!!
!! These functions are the mathematical backbone of the axisymmetric vacuum
!! response (vac_ops.calc_GH_2), so they get an external anchor:
!!
!!  - Reference values generated with mpmath 1.3 (30 significant digits),
!!    \c legenp(n-1/2,0,z,type=3) and \c legenq(n-1/2,0,z,type=3), covering
!!    the near-singular regime (z -> 1+, corresponding to nearby source and
!!    influence points in the BEM) up to the far regime (z = 100).
!!  - The same at the degrees of peeling-ballooning modes, n = 20 to 80
!!    (issue #4, catalogue item T11): Chance, Turnbull and Snyder 2006 found
!!    that toroidal functions built by upward recursion in n lose precision
!!    there; dtorh1 recurs backward for Q and forward for P and must keep
!!    its advertised precision over the same z range.
!!  - The three-term recurrence in the degree, which both P and Q satisfy:
!!    (n+1/2) F_{n+1/2} = 2 n z F_{n-1/2} - (n-1/2) F_{n-3/2},
!!    checked over a long upward sweep.
!!  - Error handling for invalid arguments.
!------------------------------------------------------------------------------!
module test_dtorh
    use testdrive, only: new_unittest, unittest_type, error_type, check, &
        &test_failed
    use num_vars, only: dp
    use str_utilities, only: i2str, r2str
    use dtorh, only: dtorh1

    implicit none
    private
    public collect_dtorh

contains
    !> Collect all tests of this suite.
    subroutine collect_dtorh(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
            &new_unittest("reference_values", test_reference_values), &
            &new_unittest("high_degree", test_high_degree), &
            &new_unittest("recurrence", test_recurrence), &
            &new_unittest("invalid_argument", test_invalid_argument) &
            &]
    end subroutine collect_dtorh

    !> Compare against high-precision reference values.
    subroutine test_reference_values(error)
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: nmax = 10                                         ! highest degree n in reference set
        ! columns: z, n, P_{n-1/2}(z), Q_{n-1/2}(z)  [mpmath, 30 digits]
        integer, parameter :: n_ref = 54
        real(dp), parameter :: ref(4,n_ref) = reshape([ &
            &1.0001_dp, 0._dp, 0.99998750035155029_dp, 6.3379714137292353_dp, &
            &1.0001_dp, 1._dp, 1.0000374994140796_dp, 4.3382633107099694_dp, &
            &1.0001_dp, 2._dp, 1.0001875041015112_dp, 3.6722723781449755_dp, &
            &1.0001_dp, 3._dp, 1.0004375369146265_dp, 3.2732653821864823_dp, &
            &1.0001_dp, 5._dp, 1.0012378519507217_dp, 2.7681187444030042_dp, &
            &1.0001_dp, 10._dp, 1.0049935972764442_dp, 2.0868932202239861_dp, &
            &1.001_dp, 0._dp, 0.99987503514404764_dp, 5.1862223889747290_dp, &
            &1.001_dp, 1._dp, 1.0003749414233338_dp, 3.1885653771920325_dp, &
            &1.001_dp, 2._dp, 1.0018754101049937_dp, 2.5269311271007230_dp, &
            &1.001_dp, 3._dp, 1.0043786919701576_dp, 2.1339936668492985_dp, &
            &1.001_dp, 5._dp, 1.0124102280785720_dp, 1.6448303650171103_dp, &
            &1.001_dp, 10._dp, 1.0504875928851103_dp, 1.0177524655942508_dp, &
            &1.01_dp, 0._dp, 0.99875350346451033_dp, 4.0316687795887199_dp, &
            &1.01_dp, 1._dp, 1.0037441576549926_dp, 2.0493184514620744_dp, &
            &1.01_dp, 2._dp, 1.0187909644872200_dp, 1.4158592547726869_dp, &
            &1.01_dp, 3._dp, 1.0441197040183519_dp, 1.0584374848354174_dp, &
            &1.01_dp, 5._dp, 1.1273059445365805_dp, 0.65142626170330841_dp, &
            &1.01_dp, 10._dp, 1.5629534659683907_dp, 0.23892979561503942_dp, &
            &1.1_dp, 0._dp, 0.98783980460580096_dp, 2.8611928721988956_dp, &
            &1.1_dp, 1._dp, 1.0369305737585080_dp, 0.97876028288694116_dp, &
            &1.1_dp, 2._dp, 1.1915515733105448_dp, 0.48178412416788182_dp, &
            &1.1_dp, 3._dp, 1.4749724247714540_dp, 0.26068388880330731_dp, &
            &1.1_dp, 5._dp, 2.6275461342503970_dp, 0.085580437597598351_dp, &
            &1.1_dp, 10._dp, 16.279079432414931_dp, 0.0067473810931886026_dp, &
            &1.5_dp, 0._dp, 0.94500633092975805_dp, 2.0189058199784232_dp, &
            &1.5_dp, 1._dp, 1.1746724294455385_dp, 0.39317514837200473_dp, &
            &1.5_dp, 2._dp, 2.0343427485811577_dp, 0.11338169008453506_dp, &
            &1.5_dp, 3._dp, 4.1776191389274553_dp, 0.036210967179681298_dp, &
            &1.5_dp, 5._dp, 21.522333339356769_dp, 0.0041745653916465122_dp, &
            &1.5_dp, 10._dp, 1836.4101325086369_dp, 2.4377561438024632e-5_dp, &
            &2.0_dp, 0._dp, 0.90128629936044730_dp, 1.6566381702365942_dp, &
            &2.0_dp, 1._dp, 1.3291381621853578_dp, 0.22401429283641564_dp, &
            &2.0_dp, 2._dp, 3.2439396660408049_dp, 0.045158724151576977_dp, &
            &2.0_dp, 3._dp, 9.5831240340193611_dp, 0.010099341583196943_dp, &
            &2.0_dp, 5._dp, 101.13072752211733_dp, 0.00057191641375056766_dp, &
            &2.0_dp, 10._dp, 50988.725501617637_dp, 5.6639453944003576e-7_dp, &
            &5.0_dp, 0._dp, 0.74574918731632961_dp, 1.0010773804561062_dp, &
            &5.0_dp, 1._dp, 2.0355638390559679_dp, 0.050629509754072014_dp, &
            &5.0_dp, 2._dp, 13.321842531267676_dp, 0.0038376048751113481_dp, &
            &5.0_dp, 3._dp, 105.35340194670783_dp, 0.00032313314844757609_dp, &
            &5.0_dp, 5._dp, 7860.4012005762194_dp, 2.5974322205521254e-6_dp, &
            &5.0_dp, 10._dp, 521569837.42284866_dp, 1.9569276952369690e-11_dp, &
            &10.0_dp, 0._dp, 0.62452096119108595_dp, 0.70380587894745619_dp, &
            &10.0_dp, 1._dp, 2.8573498347230639_dp, 0.017644903169651927_dp, &
            &10.0_dp, 2._dp, 37.889824142577157_dp, 0.00066341594620695835_dp, &
            &10.0_dp, 3._dp, 604.52277638040067_dp, 2.7713237520177744e-5_dp, &
            &10.0_dp, 5._dp, 183284.18779598399_dp, 5.4837831068257537e-8_dp, &
            &10.0_dp, 10._dp, 404445252514.10914_dp, 1.2425051713931600e-14_dp, &
            &100.0_dp, 0._dp, 0.30091748588199265_dp, 0.22214831233847302_dp, &
            &100.0_dp, 1._dp, 9.0037466610689607_dp, 0.00055538640149552274_dp, &
            &100.0_dp, 2._dp, 1200.3992489805674_dp, 2.0827532453582535e-6_dp, &
            &100.0_dp, 3._dp, 192058.47758889415_dp, 8.6783600069118692e-9_dp, &
            &100.0_dp, 5._dp, 5852908935.0997311_dp, 1.7086384492823595e-13_dp, &
            &100.0_dp, 10._dp, 1.3077930651638725e+21_dp, 3.8234264910511346e-25_dp], &
            &[4,n_ref])
        real(dp), parameter :: rtol = 1.e-11_dp                                 ! dtorh1 advertises ~1e-12 for ipre=1

        integer :: id                                                           ! counter
        integer :: n                                                            ! degree
        integer :: newn                                                         ! maximum reached degree
        integer :: ierr                                                         ! error status
        real(dp) :: z                                                           ! argument
        real(dp) :: z_prev                                                      ! previous argument
        real(dp) :: pl(0:nmax), ql(0:nmax)                                      ! toroidal harmonics

        z_prev = -1._dp
        do id = 1, n_ref
            z = ref(1,id)
            n = nint(ref(2,id))

            ! (re)calculate the whole array only when z changes (exact
            ! comparison of table values, written .gt. to satisfy
            ! -Wcompare-reals)
            if (abs(z-z_prev).gt.0._dp) then
                ierr = dtorh1(z,0,nmax,pl,ql,newn)
                call check(error, ierr, 0, 'dtorh1 failed for z = '//&
                    &trim(r2str(z)))
                if (allocated(error)) return
                call check(error, newn, nmax, 'dtorh1 did not reach nmax &
                    &for z = '//trim(r2str(z)))
                if (allocated(error)) return
                z_prev = z
            end if

            call check(error, pl(n), ref(3,id), thr=rtol, rel=.true., &
                &message='P_{n-1/2} mismatch for z = '//trim(r2str(z))//&
                &', n = '//trim(i2str(n)))
            if (allocated(error)) return
            call check(error, ql(n), ref(4,id), thr=rtol, rel=.true., &
                &message='Q_{n-1/2} mismatch for z = '//trim(r2str(z))//&
                &', n = '//trim(i2str(n)))
            if (allocated(error)) return
        end do
    end subroutine test_reference_values

    !> The same at high degree, n = 20 to 80 (T11 of issue #4).
    subroutine test_high_degree(error)
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: nmax = 80                                         ! highest degree n in reference set
        ! columns: z, n, P_{n-1/2}(z), Q_{n-1/2}(z)  [mpmath 1.3, 30 digits,
        ! legenp(n-1/2,0,z,type=3) and legenq(n-1/2,0,z,type=3)]
        integer, parameter :: n_ref = 32
        real(dp), parameter :: ref(4,n_ref) = reshape([ &
            &1.0001_dp, 20._dp, 1.020087092998752_dp, 1.4265480890409455_dp, &
            &1.0001_dp, 40._dp, 1.0815992152893484_dp, 8.2400575427848697e-1_dp, &
            &1.0001_dp, 60._dp, 1.1882453067267344_dp, 5.2537188657620365e-1_dp, &
            &1.0001_dp, 80._dp, 1.3465048568346876_dp, 3.5000539338707087e-1_dp, &
            &1.001_dp, 20._dp, 1.2100325729243437_dp, 4.9064820838809055e-1_dp, &
            &1.001_dp, 40._dp, 1.9744940464950871_dp, 1.4797664163690202e-1_dp, &
            &1.001_dp, 60._dp, 3.7904017852077835_dp, 5.0233615553527741e-2_dp, &
            &1.001_dp, 80._dp, 7.8746405869027057_dp, 1.7949252637678568e-2_dp, &
            &1.01_dp, 20._dp, 4.2383350467661752_dp, 4.2422924337499804e-2_dp, &
            &1.01_dp, 40._dp, 4.8924272869596477e+1_dp, 1.8096074674162179e-3_dp, &
            &1.01_dp, 60._dp, 6.6828392636735271e+2_dp, 8.8110674548773549e-5_dp, &
            &1.01_dp, 80._dp, 9.7281436123987432e+3_dp, 4.5360691187088299e-6_dp, &
            &1.1_dp, 20._dp, 9.5392224855221361e+2_dp, 5.7276820987365432e-5_dp, &
            &1.1_dp, 40._dp, 4.7669330585566616e+6_dp, 5.7243183450024589e-9_dp, &
            &1.1_dp, 60._dp, 2.7659455838686773e+10_dp, 6.575631366647756e-13_dp, &
            &1.1_dp, 80._dp, 1.7045101286943105e+14_dp, 8.0022337099558005e-17_dp, &
            &1.7_dp, 20._dp, 4.374561110329416e+8_dp, 4.1576388712543128e-11_dp, &
            &1.7_dp, 40._dp, 1.7578134791920325e+18_dp, 5.1727833984445275e-21_dp, &
            &1.7_dp, 60._dp, 8.1775448444531995e+27_dp, 7.4126398535815052e-31_dp, &
            &1.7_dp, 80._dp, 4.0376849078691711e+37_dp, 1.1259553790579511e-40_dp, &
            &5.0_dp, 20._dp, 3.3107299050437226e+18_dp, 1.5414033769571243e-21_dp, &
            &5.0_dp, 40._dp, 1.9047153208376743e+38_dp, 1.3396018286188375e-41_dp, &
            &5.0_dp, 60._dp, 1.2680338016351796e+58_dp, 1.3414760735283907e-61_dp, &
            &5.0_dp, 80._dp, 8.9585753936772063e+77_dp, 1.4240846228773685e-81_dp, &
            &20.0_dp, 20._dp, 2.1808988799542008e+30_dp, 5.7387640940692988e-34_dp, &
            &20.0_dp, 40._dp, 1.669274743844122e+62_dp, 3.748830699178776e-66_dp, &
            &20.0_dp, 60._dp, 1.4784113428109409e+94_dp, 2.8218703606766514e-98_dp, &
            &20.0_dp, 80._dp, 1.3895221878488605e+126_dp, 2.2517911337139819e-130_dp, &
            &100.0_dp, 20._dp, 9.4081207901064017e+43_dp, 2.6574117653097128e-48_dp, &
            &100.0_dp, 40._dp, 6.9504616014591463e+89_dp, 1.7985316460113353e-94_dp, &
            &100.0_dp, 60._dp, 5.9415312808212044e+135_dp, 1.4026266405599618e-140_dp, &
            &100.0_dp, 80._dp, 5.3899604303110644e+181_dp, 1.1596212284789459e-186_dp &
            &], [4,n_ref])
        real(dp), parameter :: rtol = 5.e-11_dp                                 ! measured worst case 1.4e-11 (Q at z = 1.0001, n = 80)

        integer :: id                                                           ! counter
        integer :: n                                                            ! degree
        integer :: newn                                                         ! maximum reached degree
        integer :: ierr                                                         ! error status
        real(dp) :: z                                                           ! argument
        real(dp) :: z_prev                                                      ! previous argument
        real(dp) :: pl(0:nmax), ql(0:nmax)                                      ! toroidal harmonics
        real(dp) :: err_P, err_Q                                                ! largest relative errors
        integer :: id_P, id_Q                                                   ! where they occur

        z_prev = -1._dp
        err_P = 0._dp
        err_Q = 0._dp
        id_P = 1
        id_Q = 1
        do id = 1, n_ref
            z = ref(1,id)
            n = nint(ref(2,id))
            if (abs(z-z_prev).gt.0._dp) then
                ierr = dtorh1(z,0,nmax,pl,ql,newn)
                call check(error, ierr, 0, 'dtorh1 failed for z = '//&
                    &trim(r2str(z)))
                if (allocated(error)) return
                call check(error, newn, nmax, 'dtorh1 did not reach nmax &
                    &for z = '//trim(r2str(z)))
                if (allocated(error)) return
                z_prev = z
            end if
            if (abs(pl(n)/ref(3,id)-1._dp).gt.err_P) then
                err_P = abs(pl(n)/ref(3,id)-1._dp)
                id_P = id
            end if
            if (abs(ql(n)/ref(4,id)-1._dp).gt.err_Q) then
                err_Q = abs(ql(n)/ref(4,id)-1._dp)
                id_Q = id
            end if
        end do
        call check(error, err_P.lt.rtol, 'P_{n-1/2}: largest relative &
            &error '//trim(r2str(err_P))//' at z = '//&
            &trim(r2str(ref(1,id_P)))//', n = '//trim(i2str(nint(ref(2,id_P)))))
        if (allocated(error)) return
        call check(error, err_Q.lt.rtol, 'Q_{n-1/2}: largest relative &
            &error '//trim(r2str(err_Q))//' at z = '//&
            &trim(r2str(ref(1,id_Q)))//', n = '//trim(i2str(nint(ref(2,id_Q)))))
    end subroutine test_high_degree

    !> Three-term recurrence in the degree:
    !! (n+1/2) F_{n+1/2} = 2 n z F_{n-1/2} - (n-1/2) F_{n-3/2}.
    subroutine test_recurrence(error)
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: nmax = 30
        real(dp), parameter :: z = 1.7_dp
        real(dp), parameter :: rtol = 1.e-10_dp

        integer :: n                                                            ! counter
        integer :: newn                                                         ! maximum reached degree
        integer :: ierr                                                         ! error status
        real(dp) :: pl(0:nmax), ql(0:nmax)                                      ! toroidal harmonics
        real(dp) :: lhs, rhs                                                    ! recurrence sides

        ierr = dtorh1(z,0,nmax,pl,ql,newn)
        call check(error, ierr, 0, 'dtorh1 failed')
        if (allocated(error)) return
        call check(error, newn, nmax, 'dtorh1 did not reach nmax')
        if (allocated(error)) return

        do n = 1, nmax-1
            ! P
            lhs = (n+0.5_dp)*pl(n+1)
            rhs = 2._dp*n*z*pl(n) - (n-0.5_dp)*pl(n-1)
            call check(error, lhs, rhs, thr=rtol, rel=.true., &
                &message='P recurrence violated at n = '//trim(i2str(n)))
            if (allocated(error)) return

            ! Q
            lhs = (n+0.5_dp)*ql(n+1)
            rhs = 2._dp*n*z*ql(n) - (n-0.5_dp)*ql(n-1)
            call check(error, lhs, rhs, thr=rtol, rel=.true., &
                &message='Q recurrence violated at n = '//trim(i2str(n)))
            if (allocated(error)) return
        end do
    end subroutine test_recurrence

    !> Arguments z <= 1 lie outside the toroidal-harmonics domain and have to
    !! be rejected.
    subroutine test_invalid_argument(error)
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: nmax = 3

        integer :: newn                                                         ! maximum reached degree
        integer :: ierr                                                         ! error status
        real(dp) :: pl(0:nmax), ql(0:nmax)                                      ! toroidal harmonics

        ierr = dtorh1(0.5_dp,0,nmax,pl,ql,newn)
        call check(error, ierr.ne.0, 'dtorh1 accepted z < 1')
        if (allocated(error)) return

        ierr = dtorh1(1.0_dp,0,nmax,pl,ql,newn)
        call check(error, ierr.ne.0, 'dtorh1 accepted z = 1')
    end subroutine test_invalid_argument
end module test_dtorh
