! (from http://patorjk.com/software/taag/ ANSI shadow)

!   ██████╗  ██████╗ ███████╗████████╗
!   ██╔══██╗██╔═══██╗██╔════╝╚══██╔══╝
!   ██████╔╝██║   ██║███████╗   ██║   
!   ██╔═══╝ ██║   ██║╚════██║   ██║   
!   ██║     ╚██████╔╝███████║   ██║   
!   ╚═╝      ╚═════╝ ╚══════╝   ╚═╝   

!------------------------------------------------------------------------------!
!   program Peeling Ballooning in 3D: postprocessing                           !
!------------------------------------------------------------------------------!
!   Author: Toon Weyens                                                        !
!   Institution: Departamento de Física,                                       !
!                Universidad Carlos III de Madrid, Spain                       !
!   Contact: tweyens@fis.uc3m.es                                               !
!------------------------------------------------------------------------------!
!   Version: 0.94                                                              !
!------------------------------------------------------------------------------!
!   References:                                                                !
!       [1] Three dimensional peeling-ballooning theory in magnetic fusion     !
!           devices, eq. (6.12) and (6.16)                                     !
!------------------------------------------------------------------------------!
#define CHCKERR if(ierr.ne.0) then; call sudden_stop(ierr); end if
program POST
    use str_ops, only: i2str
    use num_vars, only: prog_name, prog_style, ltest
    use messages, only: writo, print_goodbye, lvl_ud, print_hello, &
        &init_messages, init_time, start_time, stop_time, passed_time
    use HDF5_vars, only: init_HDF5
    use X_vars, only: init_X_vars
    use MPI_ops, only: start_MPI, stop_MPI, broadcast_input_vars
    use files_ops, only: init_files, parse_args, open_input, open_output, &
        &close_output
    use input_ops, only: read_input
    use driver_POST, only: run_driver_POST
    use test, only: generic_tests
    
    implicit none

    ! local variables
    integer :: ierr                                                             ! error
    
    !-------------------------------------------------------
    !   Initialize some routines
    !-------------------------------------------------------
    ierr = start_MPI()                                                          ! start MPI
    CHCKERR
    prog_name = 'POST'                                                          ! program name
    prog_style = 2                                                              ! post-processing part
    call print_hello                                                            ! print message with time, etc
    call init_messages                                                          ! initialize message operations
    ierr = init_files()                                                         ! initialize file operations
    CHCKERR
    call init_time                                                              ! initialize time
    call init_HDF5                                                              ! initialize HDF5
    call init_X_vars                                                            ! initialize perturbation vars
 
    !-------------------------------------------------------
    !   Read the PB3D output
    !-------------------------------------------------------
    call start_time
    call writo('Initialization')
    call lvl_ud(1)
    ierr = parse_args()                                                         ! parse argument (options are used in open_input)
    CHCKERR
    ierr = open_input()                                                         ! open the input files
    CHCKERR
    ierr = read_input()                                                         ! read input file
    CHCKERR
    ierr = open_output()                                                        ! open output file per alpha group
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
    !   Main driver
    !-------------------------------------------------------
    call start_time
    call writo('Main driver')
    call lvl_ud(1)
    ierr = run_driver_POST()
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
        use num_vars, only: rank
        use MPI_ops, only: abort_MPI
        
        ! input / output
        integer, intent(in) :: ierr                                             ! error to output
        
        ! local variables
        integer :: ierr_abort                                                   ! error to output
        
        if (ierr.ne.66) then
            call writo('>> calling routine: POST (main) of rank '//&
                &trim(i2str(rank)),persistent=.true.)
            call writo('ERROR CODE '//trim(i2str(ierr))//&
                &'. Aborting MPI rank '//trim(i2str(rank)),&
                &persistent=.true.)
            call lvl_ud(1)
            ierr_abort = abort_MPI()
        else
            ierr_abort = stop_MPI()
        end if
        if (ierr_abort.ne.0) then
            call writo('MPI cannot abort...',persistent=.true.)
            call writo('Shutting down',persistent=.true.)
        end if
        call lvl_ud(-1)
        stop
    end subroutine
end program POST
