! (from http://patorjk.com/software/taag/ ANSI shadow)

!   ██████╗ ██████╗ ██████╗ ██████╗     ██████╗  ██████╗ ███████╗████████╗
!   ██╔══██╗██╔══██╗╚════██╗██╔══██╗    ██╔══██╗██╔═══██╗██╔════╝╚══██╔══╝
!   ██████╔╝██████╔╝ █████╔╝██║  ██║    ██████╔╝██║   ██║███████╗   ██║   
!   ██╔═══╝ ██╔══██╗ ╚═══██╗██║  ██║    ██╔═══╝ ██║   ██║╚════██║   ██║   
!   ██║     ██████╔╝██████╔╝██████╔╝    ██║     ╚██████╔╝███████║   ██║   
!   ╚═╝     ╚═════╝ ╚═════╝ ╚═════╝     ╚═╝      ╚═════╝ ╚══════╝   ╚═╝   

!------------------------------------------------------------------------------!
!   Postprocessing program of Peeling Ballooning in 3D                         !
!------------------------------------------------------------------------------!
!   Author: Toon Weyens                                                        !
!   Institution: Departamento de Física,                                       !
!                Universidad Carlos III de Madrid, Spain                       !
!   Contact: tweyens@fis.uc3m.es                                               !
!------------------------------------------------------------------------------!
!   Version: 0.85                                                              !
!------------------------------------------------------------------------------!
!   References:                                                                !
!       [1] Three dimensional peeling-ballooning theory in magnetic fusion     !
!           devices, eq. (6.12) and (6.16)                                     !
!------------------------------------------------------------------------------!
#define CHCKERR if(ierr.ne.0) then; call sudden_stop(ierr); end if
program PB3D_POST
    use str_ops, only: i2str
    use num_vars, only: prog_name, prog_style, ltest, output_name
    use messages, only: writo, print_goodbye, lvl_ud, print_hello, &
        &init_messages, init_time, start_time, stop_time, passed_time
    use HDF5_ops, only: init_HDF5
    use utilities, only: init_utilities
    use MPI_ops, only: start_MPI, stop_MPI, broadcast_input_vars
    use files_ops, only: init_files, parse_args, open_input, open_output, &
        &close_output
    use input_ops, only: read_input
    use PB3D_ops, only: read_PB3D
    use driver_POST, only: run_driver_POST
    use test, only: generic_tests
    
    implicit none

    ! local variables
    integer :: ierr                                                             ! error
    write(*,*) '!!!!! FOR VMEC: COPY ALSO THE EQUILIBRIUM QUANTITIES TO BE ABLE &
        &TO CALCULATE EVERYTHING FOR A NON-ALIGNED PLOT GRID !!!'
    write(*,*) '!!!!! FOR HELENA: COPY ALSO THE FIELD-ALIGNED GRID AND &
        &INTERPOLATE THE QUANTITIES TO THIS GRID !!!'
    write(*,*) 'THIS WAY YOU HAVE TWO GRIDS: FIELD ALIGNED AND PLOT'
    
    !-------------------------------------------------------
    !   Initialize some routines
    !-------------------------------------------------------
    output_name = 'PB3D_POST_out'
    ierr = start_MPI()                                                          ! start MPI
    CHCKERR
    prog_name = 'PB3D_POST'
    prog_style = 2
    call print_hello
    call init_messages                                                          ! initialize message operations
    ierr = init_files()                                                         ! initialize file operations
    CHCKERR
    call init_utilities                                                         ! initialize utilities
    call init_time                                                              ! initialize time
    call init_HDF5                                                              ! initialize HDF5
 
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
    ierr = read_PB3D()                                                          ! read the PB3D file
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
        use num_vars, only: glb_rank
        use MPI_ops, only: abort_MPI
        
        ! input / output
        integer, intent(in) :: ierr                                             ! error to output
        
        ! local variables
        integer :: ierr_abort                                                   ! error to output
        
        if (ierr.ne.66) then
            call writo('>> calling routine: PB3D_POST (main) of rank '//&
                &trim(i2str(glb_rank)),persistent=.true.)
            call writo('ERROR CODE '//trim(i2str(ierr))//&
                &'. Aborting MPI rank '//trim(i2str(glb_rank)),&
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
end program PB3D_POST
