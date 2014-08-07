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
    use time, only: init_time, start_time, passed_time
    use test, only: test_repack, test_print_GP, test_mesh_cs, &
        &test_metric_transf, test_ang_B, test_ext_var, test_B, &
        &test_VMEC_norm_deriv, test_VMEC_conv_FHM, test_calc_RZL, &
        &test_arr_mult, test_calc_T_VF, test_calc_inv_met, test_det, &
        &test_inv, test_calc_f_deriv, test_calc_g, test_f2r, &
        &test_prepare_matrix_X, test_slepc, test_calc_interp
#if ldebug
    use num_vars, only: ltest
#endif
    use str_ops, only: r2str, i2str
    use output_ops, only: init_output_ops, lvl_ud, writo
    use VMEC_vars, only: read_VMEC                                              ! The plasma variables
    use driver, only: run_driver                                                ! Main driver
    use file_ops, only: open_input, open_output, search_file, parse_args, &
        &init_file_ops
    use input_ops, only: read_input
    use MPI_ops, only: start_MPI, stop_MPI, abort_MPI
    
    implicit none

    ! local variables
    integer :: ierr                                                             ! error
    
    !include "visitfortransimV2interface.inc"                                    ! For visit simulation support
    
    !-------------------------------------------------------
    !   Initialize some routines
    !-------------------------------------------------------
    call start_MPI(ierr)
    CHCKERR
    call init_output_ops
    call init_file_ops
    call init_time
 
    !-------------------------------------------------------
    !   Read the user-provided input file and the VMEC output
    !-------------------------------------------------------
    call start_time
    call writo('Start reading the input')
    call lvl_ud(1)
    call parse_args(ierr)                                                       ! parse arguments
    CHCKERR
    call open_input(ierr)                                                       ! open the input files
    CHCKERR
    call read_input                                                             ! read input files
    call read_VMEC(ierr)                                                        ! read VMEC file
    CHCKERR
    call open_output(ierr)                                                      ! open output file
    CHCKERR
    !call broadcast_input                                                        ! broadcast to other processors
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
        !call test_repack
        !call test_print_GP
        !call test_mesh_cs
        !call test_ext_var
        !call test_det
        !call test_inv
        !call test_f2r
        !call test_ang_B
        !call test_calc_g
        !call test_calc_f_deriv
        call test_arr_mult
        call test_calc_interp
        call test_slepc
        call test_prepare_matrix_X
        call test_metric_transf
        call test_B
        call test_calc_inv_met
        call test_calc_T_VF
        call test_calc_RZL
        call test_VMEC_conv_FHM
        call test_VMEC_norm_deriv
        call writo('')
        call passed_time
        call writo('')
        call lvl_ud(-1)
    end if
#endif

    !-------------------------------------------------------
    !   Start main driver
    !-------------------------------------------------------
    call start_time
    call writo('Start main driver')
    call lvl_ud(1)
    call run_driver(ierr)
    CHCKERR
    call writo('')
    call passed_time
    call writo('')
    call lvl_ud(-1)

    !-------------------------------------------------------
    !   cleaning up
    !-------------------------------------------------------
    call writo('Start cleaning up')
    call lvl_ud(1)
    ! CLOSE THE OUTPUT FILES, 
    call stop_MPI(ierr)
    CHCKERR
    call lvl_ud(-1)

contains
    ! stops the computations, aborting MPI, etc.
    ! as a special case, if ierr = 66, no error message is printed
    subroutine sudden_stop(ierr)
        ! input / output
        integer, intent(in) :: ierr                                             ! error to output
        
        ! local variables
        integer :: ierr_abort                                                   ! error to output
        
        if (ierr.ne.66) then
            call writo('>> calling routine: PB3D (main)')
            call writo('ERROR CODE '//trim(i2str(ierr))//&
                &'. Aborting MPI..')
            call lvl_ud(1)
            call abort_MPI(ierr_abort)
        else
            call stop_MPI(ierr_abort)
        end if
        if (ierr_abort.ne.0) then
            call writo('MPI cannot abort...')
            call writo('Shutting down')
        end if
        call lvl_ud(-1)
        stop
    end subroutine
end program PB3D 
