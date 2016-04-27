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
!                            ---                                               !
!                Department of Applied Physics                                 !
!                Eindhoven University of Technology                            !
!   Contact: weyenst@gmail.com                                                 !
!------------------------------------------------------------------------------!
!   Version: 1.18                                                              !
!------------------------------------------------------------------------------!
!   References:                                                                !
!       [1] Three dimensional peeling-ballooning theory in magnetic fusion     !
!           devices, eq. (6.12) and (6.16)                                     !
!------------------------------------------------------------------------------!
#define CHCKERR if(ierr.ne.0) then; call sudden_stop(ierr); end if
program PB3D
    use num_vars, only: prog_name, prog_style, rank, rich_restart_lvl
    use str_ops, only: r2str, i2str
    use messages, only: init_messages, lvl_ud, writo, init_time, start_time, &
        &passed_time, stop_time, print_hello, print_goodbye
    use X_vars, only: init_X_vars
    use HDF5_vars, only: init_HDF5
    use driver_eq, only: run_driver_eq
    use driver_X, only: run_driver_X
    use driver_sol, only: run_driver_sol
    use files_ops, only: open_input, open_output, parse_args, init_files, &
        &close_output
    use input_ops, only: read_input_opts, read_input_eq, print_output_in
    use input_utilities, only: dealloc_in
    use MPI_ops, only: start_MPI, stop_MPI, broadcast_input_opts, sudden_stop
    use MPI_utilities, only: wait_MPI
    use eq_ops, only: calc_normalization_const, normalize_input
    use rich_ops, only: init_rich, term_rich, do_rich, start_rich_lvl, &
        &stop_rich_lvl
    use rich_vars, only: rich_info
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
    prog_name = 'PB3D'                                                          ! program name
    prog_style = 1                                                              ! main part
    call print_hello()                                                          ! print message with time, etc
    call init_messages()                                                        ! initialize message operations
    call init_files()                                                           ! initialize file operations
    call init_time()                                                            ! initialize time
    call init_HDF5()                                                            ! initialize HDF5
    call init_X_vars()                                                          ! initialize perturbation vars
 
    !-------------------------------------------------------
    !   Read the user-provided input file and the VMEC output
    !-------------------------------------------------------
    call start_time
    call writo('Initialization')
    call lvl_ud(1)
    ierr = wait_MPI()
    CHCKERR
    if (rank.eq.0) then                                                         ! only master
        ierr = parse_args()                                                     ! parse argument (options are used in open_input)
        CHCKERR
        ierr = open_input()                                                     ! open the input files
        CHCKERR
        ierr = read_input_opts()                                                ! read input options ile
        CHCKERR
        if (rich_restart_lvl.eq.1) then
            ierr = read_input_eq()                                              ! read input equilibrium file
            CHCKERR
            call calc_normalization_const()                                     ! set up normalization constants
            call normalize_input()                                              ! normalize the input
        end if
        ierr = open_output()                                                    ! open output file
        CHCKERR
        if (rich_restart_lvl.eq.1) then
            ierr = print_output_in('in')                                        ! print input outputs
            CHCKERR
            call dealloc_in()                                                   ! clean up input from equilibrium codes
        end if
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
        ierr = generic_tests()                                                  ! run generic test
        CHCKERR
        call writo('')
        call passed_time
        call writo('')
        call lvl_ud(-1)
    end if
#endif
    
    !-------------------------------------------------------
    !   Start Richardson Extrapolation Loop
    !-------------------------------------------------------
    call start_time
    ierr = init_rich()
    CHCKERR
    call stop_time
    
    RICH: do while(do_rich())
        call start_time
        call start_rich_lvl()                                                   ! start Richardson level, setting n_r_sol and other variables
        call stop_time
        
        !-------------------------------------------------------
        !   Main Driver: Equilibrium part
        !-------------------------------------------------------
        call start_time
        call writo('Equilibrium driver'//trim(rich_info()))
        call lvl_ud(1)
        ierr = run_driver_eq()                                                  ! equilibrium driver
        CHCKERR
        call writo('')
        call passed_time
        call writo('')
        call lvl_ud(-1)
        
        !---------------------------------------------------
        !   Main driver: Perturbation part
        !---------------------------------------------------
        call start_time
        call writo('Perturbation driver'//trim(rich_info()))
        call lvl_ud(1)
        ierr = run_driver_X()                                                   ! perturbation driver
        CHCKERR
        call writo('')
        call passed_time
        call writo('')
        call lvl_ud(-1)
        
        !---------------------------------------------------
        !   Main driver: Solution part
        !---------------------------------------------------
        call start_time
        call writo('Solution driver'//trim(rich_info()))
        call lvl_ud(1)
        ierr = run_driver_sol()                                                 ! solution driver
        CHCKERR
        call writo('')
        call passed_time
        call writo('')
        call lvl_ud(-1)
        
        call start_time
        call stop_rich_lvl()
        call stop_time
    end do RICH
    
    !-------------------------------------------------------
    !   Stop Richarson Extrapolation Loop
    !-------------------------------------------------------
    call start_time
    call term_rich()                                                          ! stop Richardson loop
    call stop_time
    
    !-------------------------------------------------------
    !   Clean up
    !-------------------------------------------------------
    call writo('Clean up')
    call lvl_ud(1)
    ierr = stop_MPI()
    CHCKERR
    call close_output
    call lvl_ud(-1)
    
    call print_goodbye
end program PB3D
