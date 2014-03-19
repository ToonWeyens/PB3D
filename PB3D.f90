! (from http://patorjk.com/software/taag/)

!   ██████╗ ██████╗ ██████╗ ██████╗ 
!   ██╔══██╗██╔══██╗╚════██╗██╔══██╗
!   ██████╔╝██████╔╝ █████╔╝██║  ██║
!   ██╔═══╝ ██╔══██╗ ╚═══██╗██║  ██║
!   ██║     ██████╔╝██████╔╝██████╔╝
!   ╚═╝     ╚═════╝ ╚═════╝ ╚═════╝ 

!-------------------------------------------------------
!   Main program Peeling Ballooning in 3D
!-------------------------------------------------------
!   Author: Toon Weyens
!   Institution: Departamento de Física,
!                Universidad Carlos III de Madrid, Spain
!   Contact: tweyens@fis.uc3m.es
!-------------------------------------------------------
!   Version: 0.1
!-------------------------------------------------------
!   References: 
!       [1] Three dimensional peeling-ballooning theory in magnetic fusion
!           devices, eq. (6.12) and (6.16)
!-------------------------------------------------------
program PB3D
    use time, only: init_time, start_time, stop_time, passed_time
    use test, only: test_repack, test_write_out, test_mesh_cs, &
        &test_metric_transf, test_theta_B
    use num_vars, only: ltest
    use str_ops, only: r2str
    use output_ops, only: init_output_ops, lvl_ud, writo
    use VMEC_vars, only: read_VMEC                                              ! The plasma variables
    use driver, only: run_driver                                                ! Main driver
    use file_ops, only: open_input, open_output, search_file, parse_args, &
        &init_file_ops
    use input_ops, only: read_input

    implicit none

    !include "visitfortransimV2interface.inc"                                    ! For visit simulation support

    !-------------------------------------------------------
    !   Initialize some routines
    !-------------------------------------------------------
    call init_file_ops
    call init_output_ops
    call init_time
 
    !-------------------------------------------------------
    !   Read the user-provided input file and the VMEC output
    !-------------------------------------------------------
    call start_time
    call writo('Start reading the input')
    call lvl_ud(1)
    call parse_args                                                             ! parse arguments
    call open_input                                                             ! open the input files
    call read_input                                                             ! read input files
    call read_VMEC                                                              ! read VMEC file
    call open_output                                                            ! open output file
    call writo('')
    call passed_time
    call writo('')
    call lvl_ud(-1)

    !-------------------------------------------------------
    !   Test routines and functions
    !-------------------------------------------------------
    if (ltest) then
        call start_time
        call writo('Start tests')
        call lvl_ud(1)
        !call test_repack
        !call test_write_out
        !call test_mesh_cs
        call test_metric_transf
        call test_theta_B
        call writo('')
        call passed_time
        call writo('')
        call lvl_ud(-1)
    end if

    !-------------------------------------------------------
    !   Start main driver
    !-------------------------------------------------------
    call start_time
    call writo('Start main driver')
    call lvl_ud(1)
    call run_driver
    call writo('')
    call passed_time
    call writo('')
    call lvl_ud(-1)

    !-------------------------------------------------------
    !   cleaning up
    !-------------------------------------------------------
    call writo('Start cleaning up')
    call lvl_ud(1)
    call lvl_ud(-1)

!contains
end program PB3D 
