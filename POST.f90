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
!   Institution: ITER Organization                                             !
!   Contact: weyenst@gmail.com                                                 !
!------------------------------------------------------------------------------!
!   Version: 1.45                                                              !
!------------------------------------------------------------------------------!
!   References:                                                                !
!       [1] Three dimensional peeling-ballooning theory in magnetic fusion     !
!           devices, eq. (6.12) and (6.16)                                     !
!------------------------------------------------------------------------------!
#define CHCKERR if(ierr.ne.0) then; call sudden_stop(ierr); end if
program POST
    use str_utilities, only: i2str
    use num_vars, only: prog_name, prog_style, rank
    use messages
    use HDF5_vars, only: init_HDF5
    use MPI_ops, only: start_MPI, stop_MPI, broadcast_input_opts
    use files_ops, only: init_files, parse_args, open_input, open_output, &
        &close_output
    use input_ops, only: read_input_opts
    use driver_POST, only: run_driver_POST, init_POST
    use PB3D_ops, only: reconstruct_PB3D_in
    use input_utilities, only: dealloc_in
    use eq_utilities, only: do_eq, eq_info
    use grid_vars, only: grid_type
#if ldebug
    use num_vars, only: ltest
    use test, only: generic_tests
#endif
    
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
    call print_hello()                                                          ! print message with time, etc
    call init_output()                                                          ! initialize output
    call init_files()                                                           ! initialize file operations
    call init_time()                                                            ! initialize time
    call init_HDF5()                                                            ! initialize HDF5
 
    !-------------------------------------------------------
    !   Read the PB3D output
    !-------------------------------------------------------
    call start_time
    call writo('Initialization')
    call lvl_ud(1)
    if (rank.eq.0) then
        ierr = parse_args()                                                     ! parse argument (options are used in open_input)
        CHCKERR
        ierr = open_input()                                                     ! open the input files
        CHCKERR
        ierr = reconstruct_PB3D_in('in')                                        ! reconstruct miscellaneous PB3D output variables
        CHCKERR
        ierr = read_input_opts()                                                ! read input options file
        CHCKERR
        ierr = open_output()                                                    ! open output file per alpha group
        CHCKERR
        call dealloc_in()                                                       ! clean up input from equilibrium codes
    end if
    ierr = broadcast_input_opts()                                               ! broadcast input options to other processors
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
    call writo('Initializing POST')
    call lvl_ud(1)
    ierr = init_POST()
    CHCKERR
    call stop_time
    call writo('')
    call lvl_ud(-1)
    PAR: do while(do_eq())
        call start_time
        call writo('Main driver'//trim(eq_info()))
        call lvl_ud(1)
        ierr = run_driver_POST()
        CHCKERR
        call writo('')
        call passed_time
        call writo('')
        call lvl_ud(-1)
    end do PAR

    !-------------------------------------------------------
    !   clean up
    !-------------------------------------------------------
    call writo('Clean up')
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
