!------------------------------------------------------------------------------!
!   This module contains operations on files (open, close, etc.)               !
!------------------------------------------------------------------------------!
module file_ops
#include <PB3D_macros.h>
    use netcdf
    use num_vars, only: dp, max_str_ln, max_args
    use safe_open_mod, only: safe_open
    use str_ops, only: i2str
    use output_ops, only: lvl_ud, writo, &
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
        inc_args = [1,0,0]
    end subroutine

    ! parses the command line arguments
    ! The input arguments are saved in command_arg
    ! [MPI] only global master
    integer function parse_args() result(ierr)
        use num_vars, only: prog_name, glb_rank
        
        character(*), parameter :: rout_name = 'parse_args'
        
        ! local variables
        character(len=max_str_ln) :: open_error(3), open_help(7)                ! message for incorrect usage, help, warnings
        character(len=max_str_ln) :: open_warning(2)                            ! message for warnings
        integer :: iseq                                                         ! control the user-specified arguments
        integer :: id                                                           ! dummy integer
        
        ierr = 0                                                                ! no errors (yet)
        
        if (glb_rank.eq.0) then                                                 ! following only for master process
            ! Messages for the user
            open_error(1) = ""                                                  ! incorrect usage
            open_error(2) = "Usage: " // trim(prog_name) // &
                &" USER_INPUT VMEC_INPUT [OPTIONS]"
            open_error(3) = "Try './" // trim(prog_name) // " --help' or &
                &'./" // trim(prog_name) // " -h' for more &
                &information."
            open_help(1) = open_error(2)                                        ! help with usage
            open_help(2) = ""
            open_help(3) = "    USER_INPUT       input provided by the user"
            open_help(4) = "    VMEC_INPUT       output from VMEC, in NETCDF &
                &format"
            open_help(5) = "    [OPTIONS]        (optional) command-line &
                &options"
            open_help(6) = ""
            open_help(7) = "(both can be entered in shortened form)"
            open_warning(1) = "WARNING: More than " // trim(i2str(max_args)) &
                &// " arguments given."    ! too many input arguments
            open_warning(2) = " -> extra arguments ignored"
            
            call writo("Parsing command line arguments")
            
            ! Find number of arguments and first argument
            command_arg = ""                                                    ! has to be initialized to function correctly
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
            if (numargs.gt.max_args) then                                       ! ignore extra input arguments
                do id = 1,size(open_warning)
                    call writo(open_warning(id))
                end do
                numargs=max_args
            end if
            
            do iseq = 2, numargs
                call getcarg(iseq, command_arg(iseq), numargs)
            end do
            
            call writo("Command line arguments parsed")
        end if
    end function parse_args

    ! open the input file
    ! [MPI] Only global master
    integer function open_input() result(ierr)
        use num_vars, only: VMEC_i, input_i, ltest, glb_rank, output_name
        use VMEC_vars, only: VMEC_name
        
        character(*), parameter :: rout_name = 'open_input'
        
        ! local variables
        integer :: id                                                           ! counter
        character(len=max_str_ln) :: internal_input_error(2)                    ! internal error: returned negative file handle but not an empty name
        character(len=max_str_ln) :: internal_VMEC_error(3)                     ! internal error: returned negative file handle but not an empty name
     
        character(len=max_str_ln) :: input_exts(1)                              ! all the possible extensions of the name provided by user for the input file
        character(len=max_str_ln) :: input_con_symb(3)                          ! all the possible connectors used to check for the input file
        character(len=max_str_ln) :: VMEC_exts(1)                               ! all the possible extensions of the name provided by user for the input file
        character(len=max_str_ln) :: VMEC_con_symb(3)                           ! all the possible connectors used to check for the input file
        
        ierr = 0                                                                ! no errors (yet)
     
        if (glb_rank.eq.0) then                                                 ! following only for master process
            ! Messages for the user
            internal_input_error(1) = "input file number ok but name empty"     ! problems opening the input file
            internal_input_error(2) = "input file number negative but name&
                & not empty!"
            internal_VMEC_error(1) = "VMEC file number ok but name empty"       ! problems opening the input file
            internal_VMEC_error(2) = "VMEC file number negative but name&
                & not empty!"
            internal_VMEC_error(3) = "no VMEC file found"                       ! problems opening the input file
         
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
                    CHCKERR(internal_input_error(2))
                end if
            else
                if (input_name.eq."") then
                    ierr = 1
                    CHCKERR(internal_input_error(1))
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
                    ierr = 1
                    CHCKERR(internal_VMEC_error(3))                             ! no default for VMEC input
                else 
                    ierr = 1
                    CHCKERR(internal_VMEC_error(2))
                end if
            else
                if (VMEC_name.eq."") then
                    ierr = 1
                    CHCKERR(internal_VMEC_error(1))
                else
                    call writo('VMEC file "' // trim(VMEC_name) &
                        &// '" opened at number ' // trim(i2str(VMEC_i)))
                end if
            end if
            
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
        use output_ops, only: format_out
        
        character(*), parameter :: rout_name = 'open_output'
        
        ! local variables (also used in child functions)
        character(len=max_str_ln) :: err_msg                                    ! error message
        character(len=max_str_ln) :: full_output_name                           ! including the extension
        
        ! initialize ierr
        ierr = 0
        
        if (grp_rank.eq.0) then                                                 ! only group masters
            if (glb_rank.eq.0) call writo('Attempting to open output files')    ! but only global master outputs
            call lvl_ud(1)
            
            ! apend group number to output name if not also global master
            if (glb_rank.ne.0) output_name = trim(output_name)//'_'//&
                &trim(i2str(grp_nr))
            
            ! select output format
            select case (format_out)
                case (1)                                                        ! NETCDF
                    if (glb_rank.eq.0) call writo('Output format chosen: &
                        &NETCDF')
                    ierr = open_NETCDF()
                    CHCKERR('')
                case (2)                                                        ! matlab
                    if (glb_rank.eq.0) call writo('Output format chosen: &
                        &matlab')
                    ierr = open_matlab()
                    CHCKERR('')
                case (3)                                                        ! DISLIN
                    if (glb_rank.eq.0) call writo('Output format chosen: &
                        &DISLIN')
                    ! no need to do anything
                case (4)                                                        ! GNUplot
                    if (glb_rank.eq.0) call writo('Output format chosen: &
                        &GNUplot')
                    ierr = open_gnuplot()
                    CHCKERR('')
                case default
                    if (glb_rank.eq.0) call writo('WARNING: output format "'&
                        &// trim(i2str(format_out)) // &
                        &'" is not valid. Standard output chosen')
                    ierr = open_NETCDF()
                    CHCKERR('')
            end select
            CHCKERR('')
            call lvl_ud(-1)
            if (glb_rank.eq.0) call writo('Output files opened')
        end if
    contains
        ! Open the NETCDF file
        integer function open_NETCDF() result(ierr)
            character(*), parameter :: rout_name = 'open_NETCDF'
            
            ! initialize ierr
            ierr = 0
             
            full_output_name = trim(output_name) // '.nc'
            ierr = nf90_create(trim(full_output_name), 0, output_i)             ! 0: overwrite if exists
            if (ierr.ne.0) then
                err_msg = 'Cannot open NETCDF output file'
                CHCKERR(err_msg)
            else
                if (glb_rank.eq.0) call writo('NETCDF output file "'//&
                    &trim(full_output_name) //'" opened at number '//&
                    &trim(i2str(output_i)))
            end if
        end function open_NETCDF
            
        ! Open a .m file
        integer function open_matlab() result(ierr)
            character(*), parameter :: rout_name = 'open_matlab'
            
            ! initialize ierr
            ierr = 0
            
            output_i = n_seq_0                                                  ! start at the number indicated by n_seq_0
            full_output_name = trim(output_name) // '.m'
            call safe_open(output_i,ierr,full_output_name,'replace',&
                &'formatted',delim_in='none')
            CHCKERR('Cannot open matlab output file')
            if (glb_rank.eq.0) call writo('matlab output file "'//&
                &trim(full_output_name)//'" opened at number '//&
                &trim(i2str(output_i)))
        end function open_matlab
        
        ! Open a .dat file
        integer function open_gnuplot() result(ierr)
            character(*), parameter :: rout_name = 'open_gnuplot'
            
            ! initialize ierr
            ierr = 0
            
            output_i = n_seq_0                                                  ! start at the number indicated by n_seq_0
            full_output_name = trim(output_name) // '.dat'
            call safe_open(output_i,ierr,full_output_name,'replace',&
                &'formatted',delim_in='none')
            CHCKERR('Cannot open GNUplot output file')
            if (glb_rank.eq.0) call writo('GNUplot output file "'//&
                &trim(full_output_name)//'" opened at number '//&
                &trim(i2str(output_i)))
        end function open_gnuplot
    end function open_output
    
    ! closes the output file
    ! [MPI] only group masters
    subroutine close_output
        use num_vars, only: n_groups, grp_rank, glb_rank, output_i, &
            &output_name
        
        if (n_groups.gt.1) then
            if (glb_rank.eq.0) then
                call writo('Closing output files')
                call writo('The outputs of the other groups i are saved in &
                    &their proper output files "'//trim(output_name)//'_i"')
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
