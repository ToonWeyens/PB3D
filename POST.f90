! (from http://patorjk.com/software/taag/ ANSI shadow)

!   ██████╗  ██████╗ ███████╗████████╗
!   ██╔══██╗██╔═══██╗██╔════╝╚══██╔══╝
!   ██████╔╝██║   ██║███████╗   ██║   
!   ██╔═══╝ ██║   ██║╚════██║   ██║   
!   ██║     ╚██████╔╝███████║   ██║   
!   ╚═╝      ╚═════╝ ╚══════╝   ╚═╝   

!------------------------------------------------------------------------------!
!>  Peeling Ballooning in 3D: postprocessing
!------------------------------------------------------------------------------!
!>  \author
!!  Toon Weyens,
!!  Contact: weyenst@gmail.com
!------------------------------------------------------------------------------!
!>  \version    2.32
!!  \date       2012-2018
!!  \copyright  GNU General Public License.
!------------------------------------------------------------------------------!
!>  \see
!!  References:
!!  \cite weyens2014theory
!!  \cite Weyens2017PB3D
!------------------------------------------------------------------------------!
#define CHCKERR if(ierr.ne.0) then; call sudden_stop(ierr); end if
program POST
    use str_utilities, only: i2str
    use num_vars, only: prog_name, prog_style, rank
    use messages
    use HDF5_vars, only: init_HDF5
    use MPI_ops, only: start_MPI, stop_MPI, broadcast_input_opts, &
        &sudden_stop
    use files_ops, only: init_files, parse_args, open_input, open_output, &
        &close_output
    use input_ops, only: read_input_opts
    use driver_POST, only: run_driver_POST, init_POST, stop_POST
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
    
    !------------------------------!
    !   Initialize some routines   !
    !------------------------------!
    ierr = start_MPI()                                                          ! start MPI
    CHCKERR
    prog_name = 'POST'                                                          ! program name
    prog_style = 2                                                              ! post-processing part
    call print_hello()                                                          ! print message with time, etc
    call init_output()                                                          ! initialize output
    call init_files()                                                           ! initialize file operations
    call init_time()                                                            ! initialize time
    call init_HDF5()                                                            ! initialize HDF5
 
    !--------------------------!
    !   Read the PB3D output   !
    !--------------------------!
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
    !-------------------!
    !   Do some tests   !
    !-------------------!
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
    
    !-----------------!
    !   Main driver   !
    !-----------------!
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
    call stop_POST()

    !--------------!
    !   clean up   !
    !--------------!
    call writo('Clean up')
    call lvl_ud(1)
    ierr = stop_MPI()
    CHCKERR
    call close_output
    call lvl_ud(-1)
    
    call print_goodbye
end program POST
