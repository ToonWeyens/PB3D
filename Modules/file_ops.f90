!------------------------------------------------------------------------------!
!   This module contains operations on files (open, close, etc.)               !
!------------------------------------------------------------------------------!
module file_ops
#include <PB3D_macros.h>
    use netcdf
    use num_vars, only: dp, max_str_ln, max_args
    use safe_open_mod, only: safe_open
    use str_ops, only: i2str
    use message_ops, only: lvl_ud, writo, &
        &lvl
    implicit none
    private
    public open_input, open_output, search_file, parse_args, init_file_ops, &   ! routines
        &input_name, opt_args, close_output

    ! user-specified arguments
    integer :: numargs                                                          ! control the user-specified arguments
    character(len=max_str_ln) :: command_arg(max_args)                          ! Passed command line arguments

    ! concerning input file
    character(len=max_str_ln) :: input_name                                     ! will hold the full name of the input file

    ! options provided with command line
    character(len=max_str_ln), allocatable :: opt_args(:)
    integer, allocatable :: inc_args(:)

contains
    ! initialize the variables for the module
    ! [MPI] All ranks
    subroutine init_file_ops
        use num_vars, only: max_opts, ltest, output_name
        
        output_name = "PB3D_out"                                                ! standard output name
        ltest = .false.                                                         ! don't call the testing routines
        lvl = 1
        allocate(opt_args(max_opts), inc_args(max_opts))
        opt_args = ''
        inc_args = 0
        opt_args(1) = '-o'                                                      ! in which file to output
        opt_args(2) = '-t'
        opt_args(3) = '--test'
        opt_args(4) = '--no_guess'
        opt_args(5) = '--no_plots'
        inc_args = [1,0,0,0,0]
    end subroutine

    ! parses the command line arguments
    ! The input arguments are saved in command_arg
    ! [MPI] all processes
    integer function parse_args() result(ierr)
        use num_vars, only: prog_name
        
        character(*), parameter :: rout_name = 'parse_args'
        
        ! local variables
        character(len=max_str_ln) :: open_error(3), open_help(7)                ! message for incorrect usage, help, warnings
        character(len=max_str_ln) :: open_warning(2)                            ! message for warnings
        integer :: iseq                                                         ! control the user-specified arguments
        integer :: id                                                           ! dummy integer
        
        ierr = 0                                                                ! no errors (yet)
        
        ! Messages for the user
        open_error(1) = ""                                                      ! incorrect usage
        open_error(2) = "Usage: " // trim(prog_name) // &
            &" USER_INPUT EQUILIBRIUM_INPUT [OPTIONS]"
        open_error(3) = "Try './" // trim(prog_name) // " --help' or &
            &'./" // trim(prog_name) // " -h' for more &
            &information."
        open_help(1) = open_error(2)                                            ! help with usage
        open_help(2) = ""
        open_help(3) = "    USER_INPUT          input provided by the user"
        open_help(4) = "    EQUILIBRIUM_INPUT   output from VMEC or HELENA,&
            & in .txt format"
        open_help(5) = "    [OPTIONS]           (optional) command-line &
            &options"
        open_help(6) = ""
        open_help(7) = "(both can be entered in shortened form)"
        open_warning(1) = "WARNING: More than " // trim(i2str(max_args)) &
            &// " arguments given."    ! too many input arguments
        open_warning(2) = " -> extra arguments ignored"
        
        call writo("Parsing command line arguments")
        
        ! Find number of arguments and first argument
        command_arg = ""                                                        ! has to be initialized to function correctly
        call getcarg(1, command_arg(1), numargs)
        
        ! If help with usage requested
        if (command_arg(1).eq.'-h' .or. command_arg(1).eq.'--help') then
            do id = 1,size(open_help)
                call writo(open_help(id))
            end do
            ierr = 66
            CHCKERR('')
        endif
        
        ! If not enough arguments
        if (numargs.lt.2) then
            do id = 1,size(open_error)
                call writo(open_error(id))
            end do
            ierr = 66
            CHCKERR('')
        end if
        
        ! Get the rest of the arguments
        if (numargs.gt.max_args) then                                           ! ignore extra input arguments
            do id = 1,size(open_warning)
                call writo(open_warning(id))
            end do
            numargs=max_args
        end if
        
        do iseq = 2, numargs
            call getcarg(iseq, command_arg(iseq), numargs)
        end do
        
        call writo("Command line arguments parsed")
    end function parse_args

    ! open the input files
    ! [MPI] Only global master
    integer function open_input() result(ierr)
        use num_vars, only: eq_i, input_i, ltest, glb_rank, &
            &output_name, no_guess, no_plots, eq_style, eq_name
        
        character(*), parameter :: rout_name = 'open_input'
        
        ! local variables
        integer :: id                                                           ! counter
        character(len=max_str_ln) :: err_msg                                    ! error message
        character(len=max_str_ln) :: input_exts(1)                              ! all the possible extensions of the name provided by user for the input file
        character(len=max_str_ln) :: input_con_symb(3)                          ! all the possible connectors used to check for the input file
        character(len=max_str_ln) :: eq_exts(1)                                 ! all the possible extensions of the name provided by user for the equilibrium file
        character(len=max_str_ln) :: eq_con_symb(3)                             ! all the possible connectors used to check for the equilibrium file
        character(len=4) :: first_word                                          ! first word of equilibrium file
        
        ierr = 0                                                                ! no errors (yet)
     
        if (glb_rank.eq.0) then                                                 ! following only for master process
            call writo("Attempting to open files")
            call lvl_ud(1)
            ! check for correct input file and use default if needed
            input_name = command_arg(1)                                         ! first argument is the name of the input file
            input_exts = ["input"]
            input_con_symb = [".","-","_"]
            call search_file(input_i,input_name,input_exts,input_con_symb)
            if (input_i.lt.0) then
                if (input_name.eq."") then
                    call writo('No input file found. Default used')
                else 
                    ierr = 1
                    err_msg = 'input file number negative but name not empty!'
                    CHCKERR(err_msg)
                end if
            else
                if (input_name.eq."") then
                    ierr = 1
                    err_msg = 'input file number ok but name empty'
                    CHCKERR(err_msg)
                else
                    call writo('Input file "' // trim(input_name) &
                        &// '" opened at number ' // trim(i2str(input_i)))
                end if
            end if
            
            !  check for  equilibrium  file and  print error  if  not found  (no
            ! default!)
            eq_name = command_arg(2)
            eq_exts = ["wout"]
            eq_con_symb = [".","-","_"]
            !call search_file(eq_i,eq_name,eq_exts,eq_con_symb,'.nc')
            call search_file(eq_i,eq_name,eq_exts,eq_con_symb,'.txt')
            if (eq_i.lt.0) then
                if (eq_name.eq."") then
                    err_msg = 'no equilibrium file found'
                    ierr = 1
                    CHCKERR(err_msg)                                            ! no default for equilibrium input
                else 
                    ierr = 1
                    err_msg = 'equilibrium file number negative but name not empty!'
                    CHCKERR(err_msg)
                end if
            else
                if (eq_name.eq."") then
                    ierr = 1
                    err_msg = 'equilibrium file number ok but name empty'
                    CHCKERR(err_msg)
                else
                    call writo('equilibrium file "' // trim(eq_name) &
                        &// '" opened at number ' // trim(i2str(eq_i)))
                end if
            end if
            
            ! determine which equilibrium style (1: VMEC, 2: HELENA)
            read(eq_i,*,IOSTAT=ierr) first_word                                 ! read first word of equilibrium file
            err_msg = 'Unable to read first word of equilibrium file'
            CHCKERR(err_msg)
            if (first_word.eq.'VMEC') then                                      ! It is a VMEC file
                eq_style = 1
            else                                                                ! TEMPORARILY ASSUME IT'S A HELENA FILE
                write(*,*) '!!! BETTER IDENTIFY HELENA FILES !!!'
                eq_style = 2
            end if
            backspace(UNIT=eq_i)                                                ! go back one line in the equilibrium file
            
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
                            ierr = apply_opt(jd,id)
                            CHCKERR('')
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
        integer function apply_opt(opt_nr,arg_nr) result(ierr)
            character(*), parameter :: rout_name = 'apply_opt'
            
            ! input / output
            integer :: opt_nr, arg_nr
            
            select case(opt_nr)
                case (1)                                                        ! option -o
                    if (trim(command_arg(arg_nr+1))=='') then
                        ierr = 1
                        CHCKERR('no output file name give')
                    else
                        output_name = trim(command_arg(arg_nr+1))
                        call writo('output file "' // trim(output_name) // &
                            &'" chosen')
                    end if
                case (2,3)                                                      ! option test
#if ldebug
                    call writo('option test chosen')
                    ltest = .true.
#else
                    call writo('WARNING: option test not available. &
                        &Recompile with cpp flag ''ldebug''...')
#endif
                case (4)                                                        ! disable guessing Eigenfunction from previous Richardson level
                    call writo('option no_guess chosen: Eigenfunction not &
                        &guessed from previous Richardson level')
                    call writo('(should be used for debugging only as it &
                        &can lower performance)')
                    no_guess = .true.
                case (5)                                                        ! disable plotting
                    call writo('option no_plots chosen: plotting disabled')
                    no_plots = .true.
                case default
                    call writo('WARNING: Invalid option number')
            end select
        end function apply_opt
    end function open_input

    ! open an output file
    ! [MPI] only group masters
    integer function open_output() result(ierr)
        use num_vars, only: output_i, n_seq_0, glb_rank, grp_nr, glb_rank, &
            &grp_rank, output_name
        use message_ops, only: temp_output, temp_output_active, &
            &temp_output_id, temp_output_omitted
        
        character(*), parameter :: rout_name = 'open_output'
        
        ! local variables (also used in child functions)
        character(len=max_str_ln) :: full_output_name                           ! including the extension
        integer :: id                                                           ! counter
        
        ! initialize ierr
        ierr = 0
        
        if (grp_rank.eq.0) then                                                 ! only group masters
            if (glb_rank.eq.0) call writo('Attempting to open output files')    ! but only global master outputs
            call lvl_ud(1)
            
            ! apend group number to output name if not also global master
            if (glb_rank.ne.0) output_name = trim(output_name)//'_'//&
                &trim(i2str(grp_nr))
            
            ! open output file for the log
            output_i = n_seq_0                                                  ! start at the number indicated by n_seq_0
            full_output_name = trim(output_name) // '.txt'
            call safe_open(output_i,ierr,full_output_name,'replace',&
                &'formatted',delim_in='none')
            CHCKERR('Cannot open log output file')
            temp_output_active = .false.                                        ! no more temporary output
            
            ! write temporary output to log output file
            do id = 1,temp_output_id-1
                write(output_i,*) temp_output(id)
            end do
            
            ! print message
            if (glb_rank.eq.0) then
                call writo('log output file "'//trim(full_output_name)//&
                &'" opened at number '//trim(i2str(output_i)))
                
                call writo('temporary output copied into log file')
                
                ! check for correct completion
                if (temp_output_omitted.ne.0) call writo('WARNING: '//&
                    &trim(i2str(temp_output_omitted))//' log entries have been &
                    &omited in total')
            end if
            
            ! deallocate temporary output
            deallocate(temp_output)
            
            call lvl_ud(-1)
            if (glb_rank.eq.0) call writo('Output files opened')
        end if
    end function open_output
    
    ! closes the output file
    ! [MPI] only group masters
    subroutine close_output
        use num_vars, only: n_groups, grp_rank, glb_rank, output_i, &
            &output_name
        
        if (n_groups.gt.1) then
            if (glb_rank.eq.0) then
                call writo('The output logs of the other groups i are saved in &
                    &"'//trim(output_name)//'_i".txt')
                call writo('Closing output files')
            end if
        else
            call writo('Closing output file')
        end if
        
        if (grp_rank.eq.0) close(output_i)
    end subroutine
    
    ! looks for the full name of a file and tests for its existence
    ! output:   full name of file, empty string if non-existent
    subroutine search_file(i_unit,file_name,exts,con_symb,ext)
        use num_vars, only: n_seq_0
        
        character(len=*), intent(inout) :: file_name                            ! the name that is searched for
        character(len=*), intent(in) :: exts(:)                                 ! the possible extensions of the full_name
        character(len=*), intent(in) :: con_symb(:)                             ! the possible symbols to connect with the full_name
        character(len=*), optional, intent(in) :: ext                           ! optional extension at end of file
        integer, intent(out) :: i_unit                                          ! will hold the file handle

        character(len=max_str_ln) :: mod_file_name                               ! modified file name
        integer :: id, jd, istat

        ! try to open the given name
        i_unit = n_seq_0                                                        ! start at the number indicated by n_seq_0
        call safe_open(i_unit,istat,file_name,'old','formatted')
        if (present(ext)) mod_file_name = trim(mod_file_name) // trim(ext)      ! if an extension is provided, append it
        if (istat.eq.0) return

        ! try with the extensions provided
        call writo('Literal file "' // trim(file_name) // &
            &'" not found. Trying combinations')
        call lvl_ud(1)
        do id = 1, size(exts)                                                   ! iterate over all possible extensions
            do jd = 1, size(con_symb)                                           ! iterate over all possible connectors
                mod_file_name = trim(exts(id))//trim(con_symb(jd))&
                    &//trim(file_name) 
                if (present(ext)) mod_file_name = trim(mod_file_name) &         ! if an extension is provided, append it
                    &// trim(ext)
                call writo('trying ' // trim(mod_file_name))
                call safe_open(i_unit,istat,mod_file_name,'old','formatted')
                if (istat.eq.0) then
                    file_name = mod_file_name 
                    call lvl_ud(-1)
                    return
                end if
            end do
        end do
        ! no matches, empty file_name and i_unit = -1 returned
        file_name = ""
        i_unit = -1
        call lvl_ud(-1)
    end subroutine
end module file_ops
