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
!   Version: 0.73                                                              !
!------------------------------------------------------------------------------!
!   References:                                                                !
!       [1] Three dimensional peeling-ballooning theory in magnetic fusion     !
!           devices, eq. (6.12) and (6.16)                                     !
!------------------------------------------------------------------------------!
#define CHCKERR if(ierr.ne.0) then; call sudden_stop(ierr); end if
program PB3D
    use str_ops, only: r2str, i2str
    use messages, only: init_messages, lvl_ud, writo, init_time, &
        &start_time, passed_time, print_hello, print_goodbye
    use HDF5_ops, only: init_HDF5
    use driver, only: run_driver
    use files, only: open_input, open_output, parse_args, init_files, &
        &close_output
    use input_ops, only: read_input
    use utilities, only: init_utilities
    use MPI_ops, only: start_MPI, stop_MPI, abort_MPI, broadcast_input_vars
    use eq_ops, only: read_eq, calc_norm_const
    use test, only: generic_tests
    
    implicit none

    ! local variables
    integer :: ierr                                                             ! error
    
    write(*,*) 'CALC_V_INT SHOULD WORK WITH FAST FOURIER TRANSFORM!!!'
    
    write(*,*) 'TEST PRESSURE BALANCE'
    
    write(*,*) 'INVESTIGATE HOW IMPROVING THE DERIVATIVES CAN &
        &HELP YOU GET RID OF UNSTABLE SIDE OF SPECTRUM !!!!!!!!!!!!!!!!!!!!!!!'
    write(*,*) 'DOING THIS IN THIS ROUTINE ALREADY HELPED A LOT!!!!!'
    
    write(*,*) '!!! HDF5 PLOT CHOULD CHECK IF IT IS REALLY AXISYMMETRIC IF IT &
        &DOES A 2D PLOT BECAUSE ONE DIMENSION EQUAL TO 1 !!!'
    
    !-------------------------------------------------------
    !   Initialize some routines
    !-------------------------------------------------------
    ierr = start_MPI()                                                          ! start MPI
    CHCKERR
    call print_hello
    call init_messages                                                          ! initialize message operations
    call init_files                                                             ! initialize file operations
    call init_utilities                                                         ! initialize utilities
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
    ierr = read_input()                                                         ! read input files
    CHCKERR
    ierr = open_output()                                                        ! open output file per alpha group
    CHCKERR
    ierr = calc_norm_const()                                                    ! set up normalization constants
    CHCKERR
    ierr = broadcast_input_vars()                                               ! broadcast to other processors
    CHCKERR
    call writo('')
    call passed_time
    call writo('')
    call lvl_ud(-1)
    
    !-------------------------------------------------------
    !   Do some tests
    !-------------------------------------------------------
#if ldebug
    call start_time
    call writo('Generic Tests')
    call lvl_ud(1)
    ierr = generic_tests()
    CHCKERR
    call writo('')
    call passed_time
    call writo('')
    call lvl_ud(-1)
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
end program PB3D 
