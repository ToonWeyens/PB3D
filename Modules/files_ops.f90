!------------------------------------------------------------------------------!
!   Operations related to files                                                !
!------------------------------------------------------------------------------!
module files_ops
#include <PB3D_macros.h>
    use str_ops
    use messages
    use num_vars, only: dp, max_str_ln
    implicit none
    private
    public open_input, open_output, parse_args, init_files, &
        &input_name, opt_args, close_output

    ! user-specified arguments
    integer :: numargs                                                          ! control the user-specified arguments
    character(len=max_str_ln), allocatable :: command_arg(:)                    ! passeed command-line arguments

    ! concerning input file
    character(len=max_str_ln) :: input_name                                     ! will hold the full name of the input file

    ! options provided with command line
    character(len=max_str_ln), allocatable :: opt_args(:)
    integer, allocatable :: inc_args(:)

contains
    ! initialize the variables for the module
    ! [MPI] All ranks
    integer function init_files() result(ierr)
        use num_vars, only: ltest, output_name, prog_style
        
        character(*), parameter :: rout_name = 'init_files'
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! select according to program style
        select case (prog_style)
            case(1)                                                             ! PB3D
                output_name = "PB3D_out"                                        ! standard output name
                ltest = .false.                                                 ! don't call the testing routines
                lvl = 1
                allocate(opt_args(8), inc_args(8))
                opt_args = ''
                inc_args = 0
                opt_args(1) = '-t'
                opt_args(2) = '--test'
                opt_args(3) = '--no_guess'
                opt_args(4) = '--no_plots'
                opt_args(5) = '--no_messages'
                opt_args(6) = '-st_pc_factor_shift_type'
                opt_args(7) = '-st_pc_type'
                opt_args(8) = '-st_pc_factor_mat_solver_package'
                inc_args = [0,0,0,0,0,1,1,1]
            case(2)                                                             ! PB3D_PP
                output_name = "PB3D_out"                                        ! standard output name
                ltest = .false.                                                 ! don't call the testing routines
                lvl = 1
                allocate(opt_args(5), inc_args(5))
                opt_args = ''
                inc_args = 0
                opt_args(1) = '-t'
                opt_args(2) = '--test'
                opt_args(3) = '--no_guess'
                opt_args(4) = '--no_plots'
                opt_args(5) = '--no_messages'
                inc_args = [0,0,0,0,0]
            case default
                err_msg = 'No program style associated with '//&
                    &trim(i2str(prog_style))
                ierr = 1
                CHCKERR(err_msg)
        end select
    end function init_files

    ! parses the command line arguments
    ! The input arguments are saved in command_arg
    ! [MPI] all processes
    integer function parse_args() result(ierr)
        use num_vars, only: prog_name, prog_style
        
        character(*), parameter :: rout_name = 'parse_args'
        
        ! local variables
        character(len=max_str_ln), allocatable :: open_error(:)                 ! message for incorrect usage
        character(len=max_str_ln), allocatable :: open_help(:)                  ! message for help
        integer :: iseq                                                         ! control the user-specified arguments
        integer :: id                                                           ! dummy integer
        character(len=max_str_ln) :: dum_str                                    ! dummy string
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer :: min_args                                                     ! minimal nr. of arguments
        
        ! initialize ierr
        ierr = 0
        
        ! Messages for the user
        ! select according to program style
        select case (prog_style)
            case(1)                                                             ! PB3D
                allocate(open_error(3))
                allocate(open_help(7))
                open_error(1) = ""                                              ! incorrect usage
                open_error(2) = "Usage: " // trim(prog_name) // &
                    &" USER_INPUT EQUILIBRIUM_INPUT [OPTIONS]"
                open_error(3) = "Try './" // trim(prog_name) // " --help' or &
                    &'./" // trim(prog_name) // " -h' for more &
                    &information."
                open_help(1) = open_error(2)                                    ! help with usage
                open_help(2) = ""
                open_help(3) = "    USER_INPUT          input provided by the &
                    &user"
                open_help(4) = "    EQUILIBRIUM_INPUT   output from VMEC or &
                    &HELENA, in .txt format"
                open_help(5) = "    [OPTIONS]           (optional) &
                    &command-line options"
                open_help(6) = ""
                open_help(7) = "(both can be entered in shortened form)"
                min_args = 2
            case(2)                                                             ! PB3D_PP
                allocate(open_error(3))
                allocate(open_help(6))
                open_error(1) = ""                                              ! incorrect usage
                open_error(2) = "Usage: " // trim(prog_name) // &
                    &" PB3D_OUTPUT [OPTIONS]"
                open_error(3) = "Try './" // trim(prog_name) // " --help' or &
                    &'./" // trim(prog_name) // " -h' for more &
                    &information."
                open_help(1) = open_error(2)                                    ! help with usage
                open_help(2) = ""
                open_help(3) = "    PB3D_OUTPUT         output from PB3D code, &
                    &in .h5 format"
                open_help(4) = "    [OPTIONS]           (optional) &
                    &command-line options"
                open_help(5) = ""
                open_help(6) = "(both can be entered in shortened form)"
                min_args = 1
            case default
                err_msg = 'No program style associated with '//&
                    &trim(i2str(prog_style))
                ierr = 1
                CHCKERR(err_msg)
        end select
        
        call writo("Parsing command line arguments")
        
        ! Find number of arguments and first argument
        numargs = command_argument_count() 
        dum_str = ''
        if (numargs.gt.0) call getarg(1,dum_str)
        
        ! print error if no arguments
        if (dum_str.eq.'-h' .or. dum_str.eq.'--help') then
            do id = 1,size(open_help)
                call writo(open_help(id))
            end do
            ierr = 66                                                           ! silent stop
            CHCKERR('')
        else if (numargs.lt.min_args) then
            do id = 1,size(open_error)
                call writo(open_error(id))
            end do
            ierr = 66                                                           ! silent stop
            CHCKERR('')
        else
            allocate(command_arg(numargs))
            command_arg = ""
            command_arg(1) = dum_str
        end if
        
        ! Get the rest of the arguments
        do iseq = 2, numargs
            call getarg(iseq, command_arg(iseq))
        end do
        
        call writo("Command line arguments parsed")
    end function parse_args

    ! open the input files
    ! [MPI] Only global master
    integer function open_input() result(ierr)
        use num_vars, only: eq_i, input_i, glb_rank, prog_style, &
            &no_guess, no_plots, eq_style, eq_name, no_messages
        use files_utilities, only: search_file
#if ldebug
        use num_vars, only: ltest
#endif
        
        character(*), parameter :: rout_name = 'open_input'
        
        ! local variables
        integer :: id                                                           ! counter
        character(len=max_str_ln) :: err_msg                                    ! error message
        character(len=4) :: first_word                                          ! first word of equilibrium file (VMEC)
        integer :: first_ints(2)                                                ! first two input integers of input file (HELENA)
        integer :: istat                                                        ! status
        
        ! initialize ierr
        ierr = 0
     
        if (glb_rank.eq.0) then                                                 ! following only for master process
            call writo("Opening files")
            call lvl_ud(1)
            
            ! select depending on type
            select case (prog_style)
                case(1)                                                         ! PB3D
                    ! check for correct input file and use default if needed
                    input_name = command_arg(1)                                 ! first argument is the name of the input file
                    call search_file(input_i,input_name)
                    if (input_name.eq."") then
                        call writo('No input file found. Default used')
                    else 
                        call writo('Input file "' // trim(input_name) &
                            &// '" opened at number ' // trim(i2str(input_i)))
                    end if
                    
                    ! check for  equilibrium file and  print error if  not found
                    ! (no default!)
                    eq_name = command_arg(2)
                    call search_file(eq_i,eq_name)
                    if (eq_name.eq."") then
                        ierr = 1
                        err_msg = 'No equilibrium file found and no default &
                            &possible.'
                        CHCKERR(err_msg)
                    else
                        call writo('equilibrium file "' // trim(eq_name) &
                            &// '" opened at number ' // trim(i2str(eq_i)))
                    end if
                    
                    ! Determine which equilibrium style (1: VMEC, 2: HELENA)
                    ! set eq_style to a nonsensical value
                    eq_style = -1
                    ! Check for VMEC
                    read(eq_i,*,IOSTAT=ierr) first_word                         ! read first word of equilibrium file
                    backspace(UNIT=eq_i)                                        ! go back one line in the equilibrium file
                    err_msg = 'Unable to read first word of equilibrium file'
                    CHCKERR(err_msg)
                    if (first_word.eq.'VMEC') then                              ! It is a VMEC file
                        eq_style = 1
                    end if
                    ! Check for HELENA
                    read(eq_i,*,IOSTAT=istat) first_ints
                    backspace(UNIT=eq_i)                                        ! go back one line in the equilibrium file
                    if (istat.eq.0) then
                        eq_style = 2
                    end if
                    ! Check if something was found
                    if (eq_style.lt.1 .or. eq_style.gt.2) then
                        ierr = 1
                        call writo('Unable to recognize type of input file')
                        call writo('If you are trying to open HELENA files, &
                            &you need to have a version that outputs the &
                            &variable IAS')
                        CHCKERR('')
                    end if
                case(2)                                                         ! PB3D_PP
                    ! check  for input  file and  print error  if not  found (no
                    ! default!)
                    input_name = command_arg(1)                                 ! first argument is the name of the input file
                    call search_file(input_i,input_name)
                    if (input_name.eq."") then
                        ierr = 1
                        err_msg = 'No PB3D output file found.'
                        CHCKERR(err_msg)
                    else
                        call writo('PB3D output file "' // trim(input_name) &
                            &// '" opened at number ' // trim(i2str(eq_i)))
                    end if
                case default
                    err_msg = 'No program style associated with '//&
                        &trim(i2str(prog_style))
                    ierr = 1
                    CHCKERR(err_msg)
            end select
            
            ! set options
            if (numargs.gt.2) then                                              ! options given
                call writo('applying options')
                call lvl_ud(1)
                call read_opts
                call lvl_ud(-1)
            end if
            call lvl_ud(-1)
            call writo('Files opened')
        end if
    contains
        ! this subroutine scans for chosen options
        subroutine read_opts
            integer :: jd
            integer :: numopts
            logical :: opt_taken(size(opt_args))                                ! which of the options has been taken
            logical :: opt_found
            
            ! set number of options still to be processed and no option taken yet
            numopts = size(opt_args)
            opt_taken = .false. 
            opt_found = .false.
            
            id = 3                                                              ! start at third command argument
            do while (id.le.numargs)                                            ! try the remaining possible options
                jd = 1
                do while (jd.le.numopts .and. .not.opt_found)                   ! try all the options
                    if (trim(command_arg(id)).eq.trim(opt_args(jd))) then       ! option found
                        opt_found = .true.
                        if (opt_taken(jd)) then                                 ! option already taken
                            call writo('WARNING: Option "' // &
                                &trim(command_arg(id)) // '" already set, &
                                &ignoring...')
                        else                                                    ! option not yet taken
                            select case (prog_style)
                                case(1)                                         ! PB3D
                                    call apply_opt(jd,id)
                                case(2)                                         ! PB3D_PP
                                    call apply_opt_PP(jd)   
                                case default
                                    err_msg = 'No program style associated &
                                        &with '//trim(i2str(prog_style))
                                    ierr = 1
                                    CHCKERR(err_msg)
                            end select
                            opt_taken(jd) = .true.
                        end if
                        id = id + inc_args(jd)                                  ! skip a number of next arguments
                    end if
                    jd = jd + 1
                    ! FIND THE INCREMENT WITH BELOW FUNCTION
                end do
                if (.not.opt_found) then
                    call writo('WARNING: Option "' // &
                        &trim(command_arg(id)) // '" invalid')
                end if
                id = id + 1
                opt_found = .false.
            end do
        end subroutine
        
        ! this subroutine applies chosen options
        subroutine apply_opt(opt_nr,arg_nr)
            ! input / output
            integer :: opt_nr, arg_nr
            
            select case(opt_nr)
                case (1,2)                                                      ! option test
#if ldebug
                    call writo('option test chosen')
                    ltest = .true.
#else
                    call writo('WARNING: option test not available. &
                        &Recompile with cpp flag ''ldebug''...')
#endif
                case (3)                                                        ! disable guessing Eigenfunction from previous Richardson level
                    call writo('option no_guess chosen: Eigenfunction not &
                        &guessed from previous Richardson level')
                    no_guess = .true.
                case (4)                                                        ! disable plotting
                    call writo('option no_plots chosen: plotting disabled')
                    no_plots = .true.
                case (5)                                                        ! disable messages
                    call writo('option no_messages chosen: messages disabled')
                    no_messages = .true.
                case (6)
                    call writo('option st_pc_factor_shift_type '//&
                        &trim(command_arg(arg_nr+1))//' passed to SLEPC')
                case (7)
                    call writo('option st_pc_type '//&
                        &trim(command_arg(arg_nr+1))//' passed to SLEPC')
                case (8)
                    call writo('option st_pc_factor_mat_solver_package '//&
                        &trim(command_arg(arg_nr+1))//' passed to SLEPC')
                case default
                    call writo('WARNING: Invalid option number')
            end select
        end subroutine apply_opt
        
        subroutine apply_opt_PP(opt_nr)                                         ! postprocessing version
            ! input / output
            integer :: opt_nr
            
            select case(opt_nr)
                case (1,2)                                                      ! option test
#if ldebug
                    call writo('option test chosen')
                    ltest = .true.
#else
                    call writo('WARNING: option test not available. &
                        &Recompile with cpp flag ''ldebug''...')
#endif
                case (3)                                                        ! disable guessing Eigenfunction from previous Richardson level
                    call writo('option no_guess chosen: Eigenfunction not &
                        &guessed from previous Richardson level')
                    no_guess = .true.
                case (4)                                                        ! disable plotting
                    call writo('option no_plots chosen: plotting disabled')
                    no_plots = .true.
                case (5)                                                        ! disable messages
                    call writo('option no_messages chosen: messages disabled')
                    no_messages = .true.
                case default
                    call writo('WARNING: Invalid option number')
            end select
        end subroutine apply_opt_PP
    end function open_input

    ! open an output file
    ! [MPI] Parts by all processes, parts only by group master
    ! Note: If  this routine  is called  before the  processes are  divided into
    ! groups, the group  number is equal to the global  number, which means that
    ! only  the  global  master  will  create an  output  file,  as  should  be.
    ! Afterwards, it  can be called again  and more output files  can be created
    ! for the other groups.
    integer function open_output() result(ierr)
        use num_vars, only: output_i, grp_rank, output_name, grp_nr, n_groups
        use messages, only: temp_output, temp_output_active
        use files_utilities, only: nextunit
        use MPI_utilities, only: calc_n_groups
        
        character(*), parameter :: rout_name = 'open_output'
        
        ! local variables (also used in child functions)
        integer :: id                                                           ! counter
        character(len=max_str_ln) :: full_output_name                           ! full name
        integer :: n_groups_loc                                                 ! local n_groups
        
        ! initialize ierr
        ierr = 0
        
        ! output files common for all program styles
        if (grp_rank.eq.0) then                                                 ! only group masters
            ! user output
            call writo('Attempting to open output file')
            call lvl_ud(1)
            
            ! set local n_groups
            call calc_n_groups(n_groups_loc)
            
            ! set full output name
            full_output_name = trim(output_name)
            if (n_groups.eq.0) then                                             ! n_groups not yet initialized
                if (n_groups_loc.gt.1) full_output_name = &
                    &trim(full_output_name)//'_G1'                              ! there will be more than 1 groups: append group number (+1)
            else
                full_output_name = &
                    &trim(full_output_name)//'_G'//trim(i2str(grp_nr+1))        ! append group number (+1)
            end if
            full_output_name = trim(full_output_name)//'.txt'
            
            ! open file
            open(unit=nextunit(output_i),file=trim(full_output_name),&
                &iostat=ierr)
            CHCKERR('Cannot open log output file')
            
            ! print message
            call writo('log output file "'//trim(full_output_name)//&
                &'" opened at number '//trim(i2str(output_i)))
            
            ! if temporary output present, silently write it to log output file
            do id = 1,size(temp_output)
                write(output_i,*) temp_output(id)
            end do
            
            call lvl_ud(-1)
            call writo('Output file opened')
        else                                                                    ! not group (or at beginning global) temporary output
            ! ---------------------------------------------------------------- !
            ! Note: The temporary output of non group masters is lost...       !
            ! ---------------------------------------------------------------- !
        end if
        
        ! no more temporary output
        temp_output_active = .false.
        
        ! deallocate temporary output if allocated
        if (allocated(temp_output)) deallocate(temp_output)
    end function open_output
    
    ! closes the output file
    ! [MPI] only group masters
    subroutine close_output
        use num_vars, only: n_groups, grp_rank, output_i
        
        if (n_groups.gt.1) then
            call writo('Closing output files')
        else
            call writo('Closing output file')
        end if
        
        if (grp_rank.eq.0) close(output_i)
    end subroutine
end module files_ops
