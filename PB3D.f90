! (from http://patorjk.com/software/taag/)

!   ██████╗ ██████╗ ██████╗ ██████╗ 
!   ██╔══██╗██╔══██╗╚════██╗██╔══██╗
!   ██████╔╝██████╔╝ █████╔╝██║  ██║
!   ██╔═══╝ ██╔══██╗ ╚═══██╗██║  ██║
!   ██║     ██████╔╝██████╔╝██████╔╝
!   ╚═╝     ╚═════╝ ╚═════╝ ╚═════╝ 

!------------------------------------------------------------------------------!
!   Main program Peeling Ballooning in 3D                                      !
!------------------------------------------------------------------------------!
!   Author: Toon Weyens                                                        !
!   Institution: Departamento de Física,                                       !
!                Universidad Carlos III de Madrid, Spain                       !
!   Contact: tweyens@fis.uc3m.es                                               !
!------------------------------------------------------------------------------!
!   Version: 0.5                                                               !
!------------------------------------------------------------------------------!
!   References:                                                                !
!       [1] Three dimensional peeling-ballooning theory in magnetic fusion     !
!           devices, eq. (6.12) and (6.16)                                     !
!------------------------------------------------------------------------------!
#define CHCKERR if(ierr.ne.0) then; call sudden_stop(ierr); end if
program PB3D
    use test, only: test_repack, test_print_GP, &
        &test_metric_transf, test_ang_B, test_calc_ext_var, test_B, &
        &test_VMEC_norm_deriv, test_VMEC_conv_FHM, test_calc_RZL, &
        &test_arr_mult, test_calc_T_VF, test_calc_inv_met, test_calc_det, &
        &test_inv, test_calc_f_deriv, test_calc_g, test_fourier2real, &
        &test_prepare_X, test_slepc
#if ldebug
    use num_vars, only: ltest
#endif
    use str_ops, only: r2str, i2str
    use message_ops, only: init_message_ops, lvl_ud, writo, init_time, &
        &start_time, passed_time, print_hello, print_goodbye
    use VMEC_vars, only: read_VMEC
    use driver, only: run_driver
    use file_ops, only: open_input, open_output, search_file, parse_args, &
        &init_file_ops, close_output
    use input_ops, only: read_input
    use utilities, only: init_utilities
    use MPI_ops, only: start_MPI, stop_MPI, abort_MPI, broadcast_vars
    use eq_vars, only: calc_norm_const
    
    implicit none

    ! local variables
    integer :: ierr                                                             ! error
    
    !-------------------------------------------------------
    !   Initialize some routines
    !-------------------------------------------------------
    ierr = start_MPI()                                                          ! start MPI
    CHCKERR
    call print_hello
    call init_message_ops                                                       ! initialize message operations
    call init_file_ops                                                          ! initialize file operations
    call init_utilities                                                         ! initialize utilities
    call init_time                                                              ! initialize time
 
    !-------------------------------------------------------
    !   Read the user-provided input file and the VMEC output
    !-------------------------------------------------------
    call start_time
    call writo('Initialization')
    call lvl_ud(1)
    ierr = parse_args()                                                         ! parse argument (options are used in open_input)
    CHCKERR
    ierr = open_input()                                                         ! open the input files
    CHCKERR
    ierr = read_VMEC()                                                          ! read VMEC file
    CHCKERR
    ierr = read_input()                                                         ! read input files
    CHCKERR
    ierr = open_output()                                                        ! open output file per alpha group
    CHCKERR
    call calc_norm_const                                                        ! set up normalization constants
    ierr = broadcast_vars()                                                     ! broadcast to other processors
    CHCKERR
    call writo('')
    call passed_time
    call writo('')
    call lvl_ud(-1)

    !-------------------------------------------------------
    !   Test routines and functions
    !-------------------------------------------------------
#if ldebug
    if (ltest) then
        call start_time
        call writo('Start tests')
        call lvl_ud(1)
        ierr = test_B()
        CHCKERR
        !ierr = test_ang_B()
        !CHCKERR
        !ierr = test_repack()
        !CHCKERR
        !call test_print_GP
        !CHCKERR
        !ierr = test_calc_mesh_cs()
        !CHCKERR
        !ierr = test_calc_ext_var()
        !CHCKERR
        !ierr = test_fourier2real()
        !CHCKERR
        !ierr = test_calc_det()
        !CHCKERR
        !ierr = test_inv()
        !CHCKERR
        !ierr = test_calc_g()
        !CHCKERR
        !ierr = test_calc_f_deriv()
        !CHCKERR
        !ierr = test_arr_mult()
        !CHCKERR
        !ierr = test_metric_transf()
        !CHCKERR
        !ierr = test_slepc()
        !CHCKERR
        !ierr = test_prepare_X()
        !CHCKERR
        !ierr = test_calc_inv_met()
        !CHCKERR
        !ierr = test_calc_T_VF()
        !CHCKERR
        ierr = test_calc_RZL()
        CHCKERR
        !ierr = test_VMEC_conv_FHM()
        !CHCKERR
        !ierr = test_VMEC_norm_deriv()
        !CHCKERR
        call writo('')
        call passed_time
        call writo('')
        call lvl_ud(-1)
    end if
#endif
    
    !-------------------------------------------------------
    !   Main driver
    !-------------------------------------------------------
    call start_time
    call writo('Main driver')
    call lvl_ud(1)
    ierr = run_driver()
    CHCKERR
    call writo('')
    call passed_time
    call writo('')
    call lvl_ud(-1)

    !-------------------------------------------------------
    !   cleaning up
    !-------------------------------------------------------
    call writo('Cleanup')
    call lvl_ud(1)
    ierr = stop_MPI()
    CHCKERR
    call close_output
    call lvl_ud(-1)
    
    call print_goodbye
    
contains
    ! stops the computations, aborting MPI, etc.
    ! as a special case, if ierr = 66, no error message is printed
    subroutine sudden_stop(ierr)
        use num_vars, only: glb_rank
        
        ! input / output
        integer, intent(in) :: ierr                                             ! error to output
        
        ! local variables
        integer :: ierr_abort                                                   ! error to output
        
        if (ierr.ne.66) then
            call writo('>> calling routine: PB3D (main) of rank '//&
                &trim(i2str(glb_rank)))
            call writo('ERROR CODE '//trim(i2str(ierr))//&
                &'. Aborting MPI rank '//trim(i2str(glb_rank)))
            call lvl_ud(1)
            ierr_abort = abort_MPI()
        else
            ierr_abort = stop_MPI()
        end if
        if (ierr_abort.ne.0) then
            call writo('MPI cannot abort...')
            call writo('Shutting down')
        end if
        call lvl_ud(-1)
        stop
    end subroutine
end program PB3D 
