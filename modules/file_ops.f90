!-------------------------------------------------------
!   This module contains operations on files (open, close, etc.)
!-------------------------------------------------------
module file_ops
    use netcdf
    use num_vars, only: &
        &dp, n_seq_0, max_str_ln, max_args, max_opts, prog_name, style, max_it_r,&
        &ltest, max_it_NR, tol_NR, output_i, input_i, VMEC_i, min_alpha, &
        &max_alpha, n_alpha
    use safe_open_mod, only: safe_open
    use str_ops, only: i2str
    use output_ops, only: lvl_ud, writo, &
        &lvl, format_out
    use eq_vars, only: &
        &min_par, max_par, n_par
    use VMEC_vars, only: VMEC_name
    implicit none
    private
    public open_input, open_output, search_file, parse_args, read_input, &      ! routines
        &init_file_ops, &
        &input_name                                                             ! variables

    ! user-specified arguments
    integer :: numargs                                                          ! control the user-specified arguments
    character(len=max_str_ln) :: command_arg(max_args)                          ! Passed command line arguments

    ! concerning input file
    character(len=max_str_ln) :: input_name                                     ! will hold the full name of the input file
   
    ! concerning output file
    character(len=max_str_ln) :: output_name                                    ! will hold name of output file

    ! options provided with command line
    character(len=max_str_ln), allocatable :: opt_args(:)
    integer, allocatable :: inc_args(:)

    ! input options
    namelist /inputdata/ format_out, style, min_par, &
        &max_par, min_alpha, max_alpha, n_par, n_alpha, max_it_NR, &
        &tol_NR, max_it_r

contains
    ! initialize the variables for the module
    subroutine init_file_ops
        output_name = "PB3D_out"                                                ! standard output name
        ltest = .false.                                                         ! don't call the testing routines
        lvl = 1
        allocate(opt_args(max_opts), inc_args(max_opts))
        opt_args = ''
        inc_args = 0
        opt_args(1) = '-o'                                                      ! in which file to output
        opt_args(2) = '-t'
        opt_args(3) = '--test'
        inc_args = [1,0,0]
    end subroutine

    ! parses the command line arguments
    ! The input arguments are saved in command_arg
    subroutine parse_args()
        character(len=max_str_ln) :: open_error(3), open_help(7)                ! message for incorrect usage, help, warnings
        character(len=max_str_ln) :: open_warning(2)                            ! message for warnings
        integer :: iseq                                                         ! control the user-specified arguments
        integer :: id                                                           ! dummy integer
 
        ! Messages for the user
        open_error(1) = ""                                                      ! incorrect usage
        open_error(2) = "Usage: " // trim(prog_name) // &
            &" USER_INPUT VMEC_INPUT [OPTIONS]"
        open_error(3) = "Try './" // trim(prog_name) // " --help' for more &
            &information."
        open_help(1) = open_error(2)                                            ! help with usage
        open_help(2) = ""
        open_help(3) = "    USER_INPUT       input provided by the user"
        open_help(4) = "    VMEC_INPUT       output from VMEC, in NETCDF format"
        open_help(5) = "    [OPTIONS]        (optional) command-line options"
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
            stop
        endif

        ! If not enough arguments
        if (numargs.lt.2) then
            do id = 1,size(open_error)
                call writo(open_error(id))
            end do
            stop
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
    end subroutine

    ! open the input file
    subroutine open_input()
        integer :: id                                                           ! dummy integer
        character(len=max_str_ln) :: internal_input_error(2)                    ! internal error: returned negative file handle but not an empty name
        character(len=max_str_ln) :: internal_VMEC_error(3)                     ! internal error: returned negative file handle but not an empty name
     
        character(len=max_str_ln) :: input_exts(1)                              ! all the possible extensions of the name provided by user for the input file
        character(len=max_str_ln) :: input_con_symb(3)                          ! all the possible connectors used to check for the input file
        character(len=max_str_ln) :: VMEC_exts(1)                               ! all the possible extensions of the name provided by user for the input file
        character(len=max_str_ln) :: VMEC_con_symb(3)                           ! all the possible connectors used to check for the input file
     
        ! Messages for the user
        internal_input_error(1) = "ERROR: input file number ok but name empty!" ! problems opening the input file
        internal_input_error(2) = "ERROR: input file number negative but name&
            & not empty!"
        internal_VMEC_error(1) = "ERROR: VMEC file number ok but name empty!"   ! problems opening the input file
        internal_VMEC_error(2) = "ERROR: VMEC file number negative but name&
            & not empty!"
        internal_VMEC_error(3) = "ERROR: no VMEC file found"                    ! problems opening the input file
     
        call writo("Attempting to open files")
        call lvl_ud(1)
        ! check for correct input file and use default if needed
        input_name = command_arg(1)                                             ! first argument is the name of the input file
        input_exts = ["input"]
        input_con_symb = [".","-","_"]
        call search_file(input_i,input_name,input_exts,input_con_symb)
        if (input_i.lt.0) then
            if (input_name.eq."") then
                call writo('No input file found. Default used')
            else 
                call writo(internal_input_error(2))
                stop
            end if
        else
            if (input_name.eq."") then
                call writo(internal_input_error(1))
                stop
            else
                call writo('Input file "' // trim(input_name) &
                    &// '" opened at number ' // trim(i2str(input_i)))
            end if
        end if
        ! check for VMEC file and print error if not found (no default!)
        VMEC_name = command_arg(2)
        VMEC_exts = ["wout"]
        VMEC_con_symb = [".","-","_"]
        call search_file(VMEC_i,VMEC_name,VMEC_exts,VMEC_con_symb,'.nc')
        if (VMEC_i.lt.0) then
            if (VMEC_name.eq."") then
                call writo(internal_VMEC_error(3))                              ! no default for VMEC input
                stop
            else 
                call writo(internal_VMEC_error(2))
                stop
            end if
        else
            if (VMEC_name.eq."") then
                call writo(internal_VMEC_error(1))
                stop
            else
                call writo('VMEC file "' // trim(VMEC_name) &
                    &// '" opened at number ' // trim(i2str(VMEC_i)))
            end if
        end if
        ! set options
        if (numargs.gt.2) then                                                  ! options given
            call writo('applying options')
            call lvl_ud(1)
            call read_opts()
            call lvl_ud(-1)
        end if
        call lvl_ud(-1)
        call writo('Input files opened')
 
    contains
        ! this subroutine scans for chosen options
        subroutine read_opts()
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
                                &trim(command_arg(id)) // '" already set')
                        else                                                    ! option not yet taken
                            call apply_opt(jd,id)
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
            integer :: opt_nr, arg_nr
            
            select case(opt_nr)
                case (1)                                                        ! option -o
                    if (trim(command_arg(arg_nr+1))=='') then
                        call writo('ERROR: no output file name given')
                        stop
                    else
                        output_name = trim(command_arg(arg_nr+1))
                        call writo('output file "' // trim(output_name) // &
                            &'" chosen')
                    end if
                case (2,3)                                                        ! option test
                    call writo('option test chosen')
                    ltest = .true.
                case default
                    call writo('WARNING: Invalid option number')
            end select
        end subroutine
    end subroutine

    ! open an output file
    subroutine open_output()
        integer :: ostat
        
        call writo('Attempting to open output files')
        call lvl_ud(1)
        select case (format_out)
            case (1)                                                            ! NETCDF
                call writo('Output format chosen: NETCDF')
                call open_NETCDF
            case (2)                                                            ! matlab
                call writo('Output format chosen: matlab')
                call open_matlab
            case default
                call writo('WARNING: output format "' // &
                    &trim(i2str(format_out)) // &
                    &'" is not valid. Standard output chosen')
                call open_NETCDF
        end select
        call lvl_ud(-1)
        call writo('Output files opened')

    contains
        ! Open the NETCDF file
        subroutine open_NETCDF()
            output_name = trim(output_name) // '.nc'
            ostat = nf90_create(trim(output_name), 0, output_i)                 ! 0: overwrite if exists
            if (ostat.ne.0) then
                call writo('ERROR: Cannot open NETCDF output file')
                stop
            else
                call writo('NETCDF output file "' // trim(output_name) &
                    &//'" opened at number ' // trim(i2str(output_i)))
            end if
        end subroutine
            
        ! Open a .m file
        subroutine open_matlab()
            output_i = n_seq_0                                                    ! start at the number indicated by n_seq_0
            output_name = trim(output_name) // '.m'
            call safe_open(output_i,ostat,output_name,'replace','formatted')
            if (ostat.ne.0) then
                call writo('ERROR: Cannot open matlab output file')
                stop
            else
                call writo('matlab output file "' // trim(output_name) &
                    &//'" opened at number ' // trim(i2str(output_i)))
            end if
        end subroutine
    end subroutine

    ! reads input from user-provided input file
    subroutine read_input()
        integer :: istat                                                        ! holds status of read command
        
        if (input_i.ge.0) then                                                  ! if open_input opened a file
            call writo('Setting up user-provided input "' &       
                &// trim(input_name) // '"')
        else 
            call writo('Setting default values for input')
        end if
        call lvl_ud(1)
        
        ! initialize input variables (optionally overwritten by user later)
        if (input_i.ge.0) call writo('Initialize all the inputs with &
            &default values')
        call default_input
        
        ! read user input
        if (input_i.ge.n_seq_0) then                                            ! otherwise, defaults are loaded
            read (input_i, nml=inputdata, iostat=istat) 
            if (istat.ne.0) then
                call writo('ERROR: Error reading input file "' // &
                    &trim(input_name) // '"')
                stop
            endif
            call writo('Overwriting with user-provided file "' // &
                &trim(input_name) // '"')
        end if
        call lvl_ud(-1)
    contains
        subroutine default_input()
            use num_vars, only: pi
            
            max_it_NR = 50
            max_it_r = 5                                                           ! 3 levels of Richardson's extrapolation
            tol_NR = 1.0E-10_dp
            format_out = 1                                                      ! NETCDF output
            style = 1                                                           ! Richardson Extrapolation with normal discretization
            min_par = 0.0_dp
            max_par = 2.0_dp*pi
            n_par = 10
            min_alpha = 0.0_dp
            max_alpha = 2.0_dp*pi
            n_alpha = 10
        end subroutine
    end subroutine

    ! looks for the full name of a file and tests for its existence
    ! output:   full name of file, empty string if non-existent
    subroutine search_file(i_unit,file_name,exts,con_symb,ext)
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
