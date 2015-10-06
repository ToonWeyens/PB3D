! (from http://patorjk.com/software/taag/ ANSI shadow)

!   ██████╗ ██████╗ ██████╗ ██████╗ 
!   ██╔══██╗██╔══██╗╚════██╗██╔══██╗
!   ██████╔╝██████╔╝ █████╔╝██║  ██║
!   ██╔═══╝ ██╔══██╗ ╚═══██╗██║  ██║
!   ██║     ██████╔╝██████╔╝██████╔╝
!   ╚═╝     ╚═════╝ ╚═════╝ ╚═════╝ 

!------------------------------------------------------------------------------!
!   program Peeling Ballooning in 3D                                           !
!------------------------------------------------------------------------------!
!   Author: Toon Weyens                                                        !
!   Institution: Departamento de Física,                                       !
!                Universidad Carlos III de Madrid, Spain                       !
!   Contact: tweyens@fis.uc3m.es                                               !
!------------------------------------------------------------------------------!
!   Version: 0.92                                                              !
!------------------------------------------------------------------------------!
!   References:                                                                !
!       [1] Three dimensional peeling-ballooning theory in magnetic fusion     !
!           devices, eq. (6.12) and (6.16)                                     !
!------------------------------------------------------------------------------!
#define CHCKERR if(ierr.ne.0) then; call sudden_stop(ierr); end if
program PB3D
    use num_vars, only: ltest, prog_name, prog_style
    use str_ops, only: r2str, i2str
    use messages, only: init_messages, lvl_ud, writo, init_time, &
        &start_time, passed_time, print_hello, print_goodbye
    use HDF5_ops, only: init_HDF5
    use driver_PREP, only: run_driver_PREP
    use driver_PERT, only: run_driver_PERT
    use files_ops, only: open_input, open_output, parse_args, init_files, &
        &close_output
    use input_ops, only: read_input
    use MPI_ops, only: start_MPI, stop_MPI, broadcast_input_vars, sudden_stop
    use PB3D_ops, only: read_PB3D, reconstruct_PB3D
    use eq_ops, only: read_eq, calc_normalization_const, normalize_input
    use test, only: generic_tests
    
    implicit none

    ! local variables
    integer :: ierr                                                             ! error
    
    !-------------------------------------------------------
    !   Initialize some routines
    !-------------------------------------------------------
    ierr = start_MPI()                                                          ! start MPI
    CHCKERR
    prog_name = 'PB3D'                                                          ! program name
    prog_style = 1                                                              ! pre-perturbation part
    call print_hello                                                            ! print message with time, etc
    call init_messages                                                          ! initialize message operations
    ierr = init_files()                                                         ! initialize file operations
    CHCKERR
    call init_time                                                              ! initialize time
    call init_HDF5                                                              ! initialize HDF5
 
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
    ierr = read_eq()                                                            ! read equilibrium file
    CHCKERR
    ierr = read_input()                                                         ! read input file
    CHCKERR
    ierr = open_output()                                                        ! open output file per alpha group
    CHCKERR
    ierr = calc_normalization_const()                                           ! set up normalization constants
    CHCKERR
    ierr = normalize_input()                                                    ! normalize the input
    CHCKERR
    ierr = broadcast_input_vars()                                               ! broadcast to other processors
    CHCKERR
    call writo('')
    call passed_time
    call writo('')
    call lvl_ud(-1)
    
#if ldebug
    !-------------------------------------------------------
    !   Do some tests
    !-------------------------------------------------------
    if (ltest) then
        call start_time
        call writo('Generic Tests')
        call lvl_ud(1)
        ierr = generic_tests()
        CHCKERR
        call writo('')
        call passed_time
        call writo('')
        call lvl_ud(-1)
    end if
#endif
    
    !-------------------------------------------------------
    !   Driver: Pre-Perturbation part
    !-------------------------------------------------------
    call start_time
    call writo('Pre-perturbation driver')
    call lvl_ud(1)
    ierr = run_driver_PREP()
    CHCKERR
    call writo('')
    call passed_time
    call writo('')
    call lvl_ud(-1)
    
    !-------------------------------------------------------
    !   Read the PB3D file and broadcast its variables
    !-------------------------------------------------------
    call start_time
    call writo('Preparing perturbation driver')
    call lvl_ud(1)
    prog_style = 2                                                              ! perturbation part
    ierr = open_input()                                                         ! open the input files
    CHCKERR
    ierr = read_PB3D(.true.,.false.)                                            ! read the pre-perturbation files
    CHCKERR
    ierr = broadcast_input_vars()                                               ! broadcast to other processors
    CHCKERR
    call writo('')
    call passed_time
    call writo('')
    call lvl_ud(-1)
    
    !-------------------------------------------------------
    !   Main driver: Perturbation part
    !-------------------------------------------------------
    call start_time
    call writo('Perturbation driver')
    call lvl_ud(1)
    ierr = run_driver_PERT()
    CHCKERR
    call writo('')
    call passed_time
    call writo('')
    call lvl_ud(-1)
    
    !-------------------------------------------------------
    !   Cleaning up
    !-------------------------------------------------------
    call writo('Cleanup')
    call lvl_ud(1)
    ierr = stop_MPI()
    CHCKERR
    call close_output
    call lvl_ud(-1)
    
    call print_goodbye
end program PB3D
