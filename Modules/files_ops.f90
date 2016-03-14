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
        &opt_args, close_output, print_output_in, dealloc_in

    ! user-specified arguments
    integer :: numargs                                                          ! control the user-specified arguments
    character(len=max_str_ln), allocatable :: command_arg(:)                    ! passeed command-line arguments

    ! options provided with command line
    character(len=max_str_ln), allocatable :: opt_args(:)
    integer, allocatable :: inc_args(:)

contains
    ! initialize the variables for the module
    ! [MPI] All ranks
    integer function init_files() result(ierr)
        use num_vars, only: ltest, prog_style
        
        character(*), parameter :: rout_name = 'init_files'
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! select according to program style
        select case (prog_style)
            case(1)                                                             ! PB3D
                allocate(opt_args(13), inc_args(13))
                opt_args = ''
                inc_args = 0
                opt_args(6) = '--no_guess'
                opt_args(7) = '-st_pc_factor_shift_type'
                opt_args(8) = '-st_pc_type'
                opt_args(9) = '-st_pc_factor_mat_solver_package'
                opt_args(10) = '-eps_monitor'
                opt_args(11) = '-eps_tol'
                opt_args(12) = '-eps_ncv'
                opt_args(13) = '-eps_mpd'
                inc_args(6:13) = [0,1,1,1,0,1,1,1]
            case(2)                                                             ! POST
                allocate(opt_args(5), inc_args(5))
                opt_args = ''
                inc_args = 0
            case default
                err_msg = 'No program style associated with '//&
                    &trim(i2str(prog_style))
                ierr = 1
                CHCKERR(err_msg)
        end select
        
        ! set common option arguments
        ltest = .false.                                                         ! don't call the testing routines
        lvl = 1
        opt_args(1) = '-t'
        opt_args(2) = '--test'
        opt_args(3) = '--no_plots'
        opt_args(4) = '--no_output'
        opt_args(5) = '--no_execute_command_line'
        inc_args(1:5) = [0,0,0,0,0]
    end function init_files

    ! parses the command line arguments
    ! The input arguments are saved in command_arg
    ! [MPI] all processes
    integer function parse_args() result(ierr)
        use num_vars, only: prog_style, prog_name
        
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
                allocate(open_help(6))
                open_error(1) = ""                                              ! incorrect usage
                open_error(2) = "Usage: "//prog_name//&
                    &" USER_INPUT EQUILIBRIUM_INPUT [OPTIONS]"
                open_error(3) = "Try './"//prog_name//" --help' or './"//&
                    &prog_name//" -h' for more information."
                open_help(1) = open_error(2)                                    ! help with usage
                open_help(2) = ""
                open_help(3) = "    USER_INPUT          input provided by the &
                    &user"
                open_help(4) = "    EQUILIBRIUM_INPUT   output from VMEC or &
                    &HELENA, in .txt format"
                open_help(5) = "    [OPTIONS]           (optional) &
                    &command-line options"
                open_help(6) = ""
                min_args = 2
            case(2)                                                             ! POST
                allocate(open_error(3))
                allocate(open_help(6))
                open_error(1) = ""                                              ! incorrect usage
                open_error(2) = "Usage: "//prog_name//&
                    &" USER_INPUT PB3D_OUTPUT [OPTIONS]"
                open_error(3) = "Try './"//prog_name//" --help' or './"//&
                    &prog_name//" -h' for more information."
                open_help(1) = open_error(2)                                    ! help with usage
                open_help(2) = ""
                open_help(3) = "    USER_INPUT          input provided by the &
                    &user"
                open_help(4) = "    PB3D_OUTPUT    output from PB3D code, &
                    &in .h5 format"
                open_help(5) = "    [OPTIONS]           (optional) &
                    &command-line options"
                open_help(6) = ""
                min_args = 2
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
    ! [MPI] only master
    integer function open_input() result(ierr)
        use num_vars, only: eq_i, input_i, rank, prog_style, no_plots, &
            &eq_style, eq_name, no_output, PB3D_i, PB3D_name, input_name, &
            &no_execute_command_line, output_name, prog_name
        use files_utilities, only: search_file
        use rich_vars, only: no_guess
#if ldebug
        use num_vars, only: ltest
#endif
        
        character(*), parameter :: rout_name = 'open_input'
        
        ! local variables
        integer :: id                                                           ! counter
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer :: first_ints(2)                                                ! first two input integers of input file (HELENA)
        integer :: istat                                                        ! status
        character(len=max_str_ln) :: file_ext                                   ! file extension
        
        ! initialize ierr
        ierr = 0
     
        if (rank.eq.0) then                                                     ! following only for master process
            call writo("Opening files")
            call lvl_ud(1)
            
            ! select depending on program style
            select case (prog_style)
                case(1)                                                         ! PB3D
                    ! check for correct input file and use default if needed
                    input_name = command_arg(1)                                 ! first argument is the name of the input file
                    call search_file(input_i,input_name)
                    if (input_name.eq."") then
                        call writo('No input file found. Default used')
                    else 
                        call writo('Input file "' // trim(input_name) &
                            &// '" opened')
                    end if
                    
                    ! check for  equilibrium file and  print error if  not found
                    ! (no default!)
                    eq_name = command_arg(2)
                    call search_file(eq_i,eq_name)
                    if (eq_name.eq."") then
                        ierr = 1
                        err_msg = 'No equilibrium file found.'
                        CHCKERR(err_msg)
                    else
                        call writo('equilibrium file "' // trim(eq_name) &
                            &// '" opened')
                    end if
                    
                    ! Determine which equilibrium style (1: VMEC, 2: HELENA)
                    ! set eq_style to a nonsensical value
                    eq_style = -1
                    
                    ! Check for VMEC
                    do id = len(eq_name)-1,1,-1
                        if (eq_name(id:id).eq.'.') then
                            file_ext = eq_name(id+1:len(eq_name))
                            exit
                        end if
                    end do
                    if (trim(file_ext).eq.'nc') then
                        eq_style = 1
                        close(eq_i)
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
                    
                    ! set PB3D output name
                    PB3D_name = prog_name//'_'//trim(output_name)//'.h5'
                case(2)                                                         ! POST
                    ! check for correct input file and use default if needed
                    input_name = command_arg(1)                                 ! first argument is the name of the input file
                    call search_file(input_i,input_name)
                    if (input_name.eq."") then
                        call writo('No input file found. Default used')
                    else 
                        call writo('Input file "' // trim(input_name) &
                            &// '" opened')
                    end if
                    
                    ! check for PB3D_PB3D file and  print error if not found (no
                    ! default!)
                    PB3D_name = command_arg(2)
                    call search_file(PB3D_i,PB3D_name)
                    if (PB3D_name.eq."") then
                        ierr = 1
                        err_msg = 'No PB3D file found.'
                        CHCKERR(err_msg)
                    else
                        call writo('PB3D output file "' // trim(PB3D_name) &
                            &// '" opened')
                        close(PB3D_i)
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
                            select case(jd)
                                ! common options 1..5
                                case (1,2)                                      ! option test
#if ldebug
                                    call writo('option test chosen')
                                    ltest = .true.
#else
                                    call writo('WARNING: option test not &
                                        &available. Recompile with cpp &
                                        &flag ''ldebug''...')
#endif
                                case (3)                                        ! disable plotting
                                    call writo('option no_plots chosen: &
                                        &plotting disabled')
                                    no_plots = .true.
                                case (4)                                        ! disable messages
                                    call writo('option no_output chosen: &
                                        &messages disabled')
                                    no_output = .true.
                                case (5)                                        ! disable execute_command_line
                                    call writo('option no_execute_command_line &
                                        &chosen: execute_command_line disabled')
                                    no_execute_command_line = .true.
                                ! specific options for each program style
                                case default
                                    select case (prog_style)
                                        case(1)                                 ! PB3D
                                            call apply_opt_PB3D(jd,id)
                                        case(2)                                 ! POST
                                            call writo('WARNING: Invalid &
                                                &option number')
                                        case default
                                            err_msg = 'No program style &
                                                &associated with '//&
                                                &trim(i2str(prog_style))
                                            ierr = 1
                                            CHCKERR(err_msg)
                                    end select
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
        subroutine apply_opt_PB3D(opt_nr,arg_nr)                                ! PB3D version
            ! input / output
            integer :: opt_nr, arg_nr
            
            select case(opt_nr)
                case (6)                                                        ! disable guessing Eigenfunction from previous Richardson level
                    call writo('option no_guess chosen: Eigenfunction not &
                        &guessed from previous Richardson level')
                    no_guess = .true.
                case (7)
                    call writo('option st_pc_factor_shift_type '//&
                        &trim(command_arg(arg_nr+1))//' passed to SLEPC')
                case (8)
                    call writo('option st_pc_type '//&
                        &trim(command_arg(arg_nr+1))//' passed to SLEPC')
                case (9)
                    call writo('option st_pc_factor_mat_solver_package '//&
                        &trim(command_arg(arg_nr+1))//' passed to SLEPC')
                case (10)
                    call writo('option eps_monitor passed to SLEPC')
                case (11)
                    call writo('option eps_tol passed to SLEPC')
                case (12)
                    call writo('option eps_ncv passed to SLEPC')
                case (13)
                    call writo('option eps_mpd passed to SLEPC')
                case default
                    call writo('WARNING: Invalid option number')
            end select
        end subroutine apply_opt_PB3D
    end function open_input

    ! Open an output file and write (PB3D) or read (POST) the common variables.
    ! Also sets some output variables.
    ! Note: If  there is Richardson restart,  no HDF5 files are  opened for PB3D
    ! and the output log file name is different.
    integer function open_output() result(ierr)
        use num_vars, only: prog_style, output_i, output_name, prog_name, &
            &rich_restart_lvl
        use messages, only: temp_output, temp_output_active
        use files_utilities, only: nextunit
        use HDF5_ops, only: create_output_HDF5, print_HDF5_arrs
        use HDF5_vars, only: var_1D_type
        
        character(*), parameter :: rout_name = 'open_output'
        
        ! local variables (also used in child functions)
        integer :: id                                                           ! counter
        character(len=max_str_ln) :: full_output_name                           ! full name
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! user output
        call writo('Attempting to open output files')
        call lvl_ud(1)
        
        ! append extension to output name
        full_output_name = prog_name//'_'//trim(output_name)//'.txt'
        
        ! actions depending on Richardson restart level
        if (rich_restart_lvl.eq.0) then
            ! open file after wiping it
            open(unit=nextunit(output_i),file=trim(full_output_name),&
                &status='replace',iostat=ierr)
            
            ! print message
            call writo('log output file "'//trim(full_output_name)//&
                &'" opened at number '//trim(i2str(output_i)))
        else
            ! append to existing file
            open(unit=nextunit(output_i),file=trim(full_output_name),&
                &status='old',position='append',iostat=ierr)
            
            ! print message
            call writo('log output file "'//trim(full_output_name)//&
                &'" reopened at number '//trim(i2str(output_i)))
        end if
        
        ! if temporary output present, silently write it to log output file
        do id = 1,size(temp_output)
            write(output_i,*) temp_output(id)
        end do
        
        ! specific actions for program styles
        select case (prog_style)
            case (1)                                                            ! PB3D
                ! create HDF5 file for output if no restart
                if (rich_restart_lvl.eq.0) then
                    ierr = create_output_HDF5()
                    CHCKERR('')
                end if
            case (2)                                                            ! POST
                ! do nothing
            case default
                ierr = 1
                err_msg = 'No program style associated with '//&
                    &trim(i2str(prog_style))
                CHCKERR(err_msg)
        end select
        
        ! no more temporary output
        temp_output_active = .false.
        
        ! deallocate temporary output if allocated
        if (allocated(temp_output)) deallocate(temp_output)
        
        call lvl_ud(-1)
        call writo('Output files opened')
    end function open_output
    
    ! Print input quantities to an output file:
    !   - misc_in:   prog_version, eq_style, rho_style, R_0, pres_0, B_0, psi_0
    !                rho_0, T_0, vac_perm, use_pol_flux_F, use_pol_flux_E,
    !                use_normalization, norm_disc_prec_eq, n_r_in, n_par_X
    !   - misc_in_V: is_asym_V, is_freeb_V, mnmax_V, mpol_V, ntor_V, gam_V
    !   - flux_q_H:  flux_t_V, Dflux_t_V, flux_p_V, Dflux_p_V, pres_V, rot_t_V
    !   - mn_V
    !   - RZL_V:     R_V_c, R_V_s, Z_V_c, Z_V_s, L_V_c, L_V_s
    !   - B_V_sub:   B_V_sub_c, B_V_sub_s, B_V_c, B_V_s, jac_V_c, jac_V_s
    !   - misc_in_H: ias, nchi
    !   - RZ_H:      R_H, Z_H
    !   - chi_H
    !   - qs_H
    !   - RBphi_H
    !   - pres_H
    !   - flux_p_H
    !   - misc_X:    prim_X, n_mod_X, min_sec_X, max_sec_X, norm_disc_prec_X,
    !                norm_style, U_style, X_style, matrix_SLEPC_style
    !   - misc_sol:  min_r_sol, max_r_sol, alpha, norm_disc_prec_sol, BC_style,
    !                EV_BC, EV_BC
    integer function print_output_in() result(ierr)
        use num_vars, only: eq_style, rho_style, prog_version, use_pol_flux_E, &
            &use_pol_flux_F, use_normalization, norm_disc_prec_eq, PB3D_name, &
            &norm_disc_prec_X, norm_style, U_style, X_style, tol_norm, &
            &matrix_SLEPC_style, BC_style, EV_style, norm_disc_prec_sol, EV_BC
        use eq_vars, only: R_0, pres_0, B_0, psi_0, rho_0, T_0, vac_perm, &
            &max_flux_E, max_flux_F
        use grid_vars, onLy: n_r_in, n_r_eq, n_par_X
        use grid_ops, only: calc_norm_range
        use X_vars, only: min_r_sol, max_r_sol, min_sec_X, max_sec_X, prim_X, &
            &n_mod_X
        use sol_vars, only: alpha
        use HDF5_ops, only: print_HDF5_arrs
        use HDF5_vars, only: var_1D_type, &
            &max_dim_var_1D
        use HELENA_vars, only: chi_H, flux_p_H, R_H, Z_H, nchi, ias, qs_H, &
            &pres_H, RBphi_H
        use VMEC, only: is_freeb_V, mnmax_V, mpol_V, ntor_V, is_asym_V, gam_V, &
            &R_V_c, R_V_s, Z_V_c, Z_V_s, L_V_c, L_V_s, mnmax_V, mn_V, rot_t_V, &
            &pres_V, flux_t_V, Dflux_t_V, flux_p_V, Dflux_p_V, nfp_V
#if ldebug
        use HELENA_vars, only: h_H_11, h_H_12, h_H_33
        use VMEC, only: B_V_sub_c, B_V_sub_s, B_V_c, B_V_s, jac_V_c, jac_V_s
#endif
        
        character(*), parameter :: rout_name = 'print_output_in'
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        type(var_1D_type), allocatable, target :: in_1D(:)                      ! 1D equivalent of input variables
        type(var_1D_type), pointer :: in_1D_loc => null()                       ! local element in in_1D
        integer :: id                                                           ! counter
        integer :: in_limits(2)                                                 ! min. and max. index of input variable grid of this process
        
        ! initialize ierr
        ierr = 0
        
        ! user output
        call writo('Writing input variables to output file')
        call lvl_ud(1)
        
        ! user output
        call writo('Preparing variables for writing')
        call lvl_ud(1)
        
        ! calculate limits of input range
        ierr = calc_norm_range(in_limits=in_limits)
        CHCKERR('')
        n_r_eq = in_limits(2)-in_limits(1)+1
        call writo('Only the normal range '//trim(i2str(in_limits(1)))//'..'//&
            &trim(i2str(in_limits(2)))//' is written')
        call writo('(input range relevant to the solution range '//&
            &trim(r2strt(min_r_sol))//'..'//trim(r2strt(max_r_sol))//&
            &' with tolerance '//trim(r2strt(tol_norm))//')')
        
        ! set up 1D equivalents of input variables
        allocate(in_1D(max_dim_var_1D))
        
        ! Set up common variables in_1D
        id = 1
        
        ! misc_in
        in_1D_loc => in_1D(id); id = id+1
        in_1D_loc%var_name = 'misc_in'
        allocate(in_1D_loc%tot_i_min(1),in_1D_loc%tot_i_max(1))
        allocate(in_1D_loc%loc_i_min(1),in_1D_loc%loc_i_max(1))
        in_1D_loc%loc_i_min = [1]
        in_1D_loc%loc_i_max = [19]
        in_1D_loc%tot_i_min = in_1D_loc%loc_i_min
        in_1D_loc%tot_i_max = in_1D_loc%loc_i_max
        allocate(in_1D_loc%p(19))
        in_1D_loc%p = [prog_version,eq_style*1._dp,rho_style*1._dp,&
            &R_0,pres_0,B_0,psi_0,rho_0,T_0,vac_perm,-1._dp,-1._dp,&
            &-1._dp,norm_disc_prec_eq*1._dp,n_r_in*1._dp,n_r_eq*1._dp,&
            &n_par_X*1._dp,max_flux_E,max_flux_F]
        if (use_pol_flux_E) in_1D_loc%p(11) = 1._dp
        if (use_pol_flux_F) in_1D_loc%p(12) = 1._dp
        if (use_normalization) in_1D_loc%p(13) = 1._dp
        
        ! misc_in_V or misc_in_H, depending on equilibrium style
        select case (eq_style)
            case (1)                                                            ! VMEC
                ! misc_in_V
                in_1D_loc => in_1D(id); id = id+1
                in_1D_loc%var_name = 'misc_in_V'
                allocate(in_1D_loc%tot_i_min(1),&
                    &in_1D_loc%tot_i_max(1))
                allocate(in_1D_loc%loc_i_min(1),&
                    &in_1D_loc%loc_i_max(1))
                in_1D_loc%loc_i_min = [1]
                in_1D_loc%loc_i_max = [7]
                in_1D_loc%tot_i_min = in_1D_loc%loc_i_min
                in_1D_loc%tot_i_max = in_1D_loc%loc_i_max
                allocate(in_1D_loc%p(7))
                in_1D_loc%p = [-1._dp,-1._dp,mnmax_V*1._dp,mpol_V*1._dp,&
                    &ntor_V*1._dp,nfp_V*1._dp,gam_V]
                if (is_asym_V) in_1D_loc%p(1) = 1._dp
                if (is_freeb_V) in_1D_loc%p(2) = 1._dp
                
                ! flux_t_V
                in_1D_loc => in_1D(id); id = id+1
                in_1D_loc%var_name = 'flux_t_V'
                allocate(in_1D_loc%tot_i_min(1),in_1D_loc%tot_i_max(1))
                allocate(in_1D_loc%loc_i_min(1),in_1D_loc%loc_i_max(1))
                in_1D_loc%loc_i_min = [1]
                in_1D_loc%loc_i_max = [n_r_eq]
                in_1D_loc%tot_i_min = in_1D_loc%loc_i_min
                in_1D_loc%tot_i_max = in_1D_loc%loc_i_max
                allocate(in_1D_loc%p(n_r_eq))
                in_1D_loc%p = flux_t_V(in_limits(1):in_limits(2))
                
                ! Dflux_t_V
                in_1D_loc => in_1D(id); id = id+1
                in_1D_loc%var_name = 'Dflux_t_V'
                allocate(in_1D_loc%tot_i_min(1),in_1D_loc%tot_i_max(1))
                allocate(in_1D_loc%loc_i_min(1),in_1D_loc%loc_i_max(1))
                in_1D_loc%loc_i_min = [1]
                in_1D_loc%loc_i_max = [n_r_eq]
                in_1D_loc%tot_i_min = in_1D_loc%loc_i_min
                in_1D_loc%tot_i_max = in_1D_loc%loc_i_max
                allocate(in_1D_loc%p(n_r_eq))
                in_1D_loc%p = Dflux_t_V(in_limits(1):in_limits(2))
                
                ! flux_p_V
                in_1D_loc => in_1D(id); id = id+1
                in_1D_loc%var_name = 'flux_p_V'
                allocate(in_1D_loc%tot_i_min(1),in_1D_loc%tot_i_max(1))
                allocate(in_1D_loc%loc_i_min(1),in_1D_loc%loc_i_max(1))
                in_1D_loc%loc_i_min = [1]
                in_1D_loc%loc_i_max = [n_r_eq]
                in_1D_loc%tot_i_min = in_1D_loc%loc_i_min
                in_1D_loc%tot_i_max = in_1D_loc%loc_i_max
                allocate(in_1D_loc%p(n_r_eq))
                in_1D_loc%p = flux_p_V(in_limits(1):in_limits(2))
                
                ! Dflux_p_V
                in_1D_loc => in_1D(id); id = id+1
                in_1D_loc%var_name = 'Dflux_p_V'
                allocate(in_1D_loc%tot_i_min(1),in_1D_loc%tot_i_max(1))
                allocate(in_1D_loc%loc_i_min(1),in_1D_loc%loc_i_max(1))
                in_1D_loc%loc_i_min = [1]
                in_1D_loc%loc_i_max = [n_r_eq]
                in_1D_loc%tot_i_min = in_1D_loc%loc_i_min
                in_1D_loc%tot_i_max = in_1D_loc%loc_i_max
                allocate(in_1D_loc%p(n_r_eq))
                in_1D_loc%p = Dflux_p_V(in_limits(1):in_limits(2))
                
                ! pres_V
                in_1D_loc => in_1D(id); id = id+1
                in_1D_loc%var_name = 'pres_V'
                allocate(in_1D_loc%tot_i_min(1),in_1D_loc%tot_i_max(1))
                allocate(in_1D_loc%loc_i_min(1),in_1D_loc%loc_i_max(1))
                in_1D_loc%loc_i_min = [1]
                in_1D_loc%loc_i_max = [n_r_eq]
                in_1D_loc%tot_i_min = in_1D_loc%loc_i_min
                in_1D_loc%tot_i_max = in_1D_loc%loc_i_max
                allocate(in_1D_loc%p(n_r_eq))
                in_1D_loc%p = pres_V(in_limits(1):in_limits(2))
                
                ! rot_t_V
                in_1D_loc => in_1D(id); id = id+1
                in_1D_loc%var_name = 'rot_t_V'
                allocate(in_1D_loc%tot_i_min(1),in_1D_loc%tot_i_max(1))
                allocate(in_1D_loc%loc_i_min(1),in_1D_loc%loc_i_max(1))
                in_1D_loc%loc_i_min = [1]
                in_1D_loc%loc_i_max = [n_r_eq]
                in_1D_loc%tot_i_min = in_1D_loc%loc_i_min
                in_1D_loc%tot_i_max = in_1D_loc%loc_i_max
                allocate(in_1D_loc%p(n_r_eq))
                in_1D_loc%p = rot_t_V(in_limits(1):in_limits(2))
                
                ! mn_V
                in_1D_loc => in_1D(id); id = id+1
                in_1D_loc%var_name = 'mn_V'
                allocate(in_1D_loc%tot_i_min(2),in_1D_loc%tot_i_max(2))
                allocate(in_1D_loc%loc_i_min(2),in_1D_loc%loc_i_max(2))
                in_1D_loc%loc_i_min = [1,1]
                in_1D_loc%loc_i_max = [mnmax_V,2]
                in_1D_loc%tot_i_min = in_1D_loc%loc_i_min
                in_1D_loc%tot_i_max = in_1D_loc%loc_i_max
                allocate(in_1D_loc%p(size(mn_V)))
                in_1D_loc%p = reshape(mn_V,[size(mn_V)])
                
                ! RZL_V
                in_1D_loc => in_1D(id); id = id+1
                in_1D_loc%var_name = 'RZL_V'
                allocate(in_1D_loc%tot_i_min(3),in_1D_loc%tot_i_max(3))
                allocate(in_1D_loc%loc_i_min(3),in_1D_loc%loc_i_max(3))
                in_1D_loc%loc_i_min = [1,1,1]
                in_1D_loc%loc_i_max = [mnmax_V,n_r_eq,6]
                in_1D_loc%tot_i_min = in_1D_loc%loc_i_min
                in_1D_loc%tot_i_max = in_1D_loc%loc_i_max
                allocate(in_1D_loc%p(6*mnmax_V*n_r_eq))
                in_1D_loc%p = reshape([R_V_c(:,in_limits(1):in_limits(2),0),&
                    &R_V_s(:,in_limits(1):in_limits(2),0),&
                    &Z_V_c(:,in_limits(1):in_limits(2),0),&
                    &Z_V_s(:,in_limits(1):in_limits(2),0),&
                    &L_V_c(:,in_limits(1):in_limits(2),0),&
                    &L_V_s(:,in_limits(1):in_limits(2),0)],&
                    &[6*mnmax_V*n_r_eq]) 
                
#if ldebug
                ! B_V_sub
                in_1D_loc => in_1D(id); id = id+1
                in_1D_loc%var_name = 'B_V_sub'
                allocate(in_1D_loc%tot_i_min(4),in_1D_loc%tot_i_max(4))
                allocate(in_1D_loc%loc_i_min(4),in_1D_loc%loc_i_max(4))
                in_1D_loc%loc_i_min = [1,1,1,1]
                in_1D_loc%loc_i_max = [mnmax_V,n_r_eq,3,2]
                in_1D_loc%tot_i_min = in_1D_loc%loc_i_min
                in_1D_loc%tot_i_max = in_1D_loc%loc_i_max
                allocate(in_1D_loc%p(2*mnmax_V*n_r_eq*3))
                in_1D_loc%p = reshape([&
                    &B_V_sub_c(:,in_limits(1):in_limits(2),:),&
                    &B_V_sub_s(:,in_limits(1):in_limits(2),:)],&
                    &[2*mnmax_V*n_r_eq*3])
                
                ! B_V
                in_1D_loc => in_1D(id); id = id+1
                in_1D_loc%var_name = 'B_V'
                allocate(in_1D_loc%tot_i_min(3),in_1D_loc%tot_i_max(3))
                allocate(in_1D_loc%loc_i_min(3),in_1D_loc%loc_i_max(3))
                in_1D_loc%loc_i_min = [1,1,1]
                in_1D_loc%loc_i_max = [mnmax_V,n_r_eq,2]
                in_1D_loc%tot_i_min = in_1D_loc%loc_i_min
                in_1D_loc%tot_i_max = in_1D_loc%loc_i_max
                allocate(in_1D_loc%p(2*mnmax_V*n_r_eq))
                in_1D_loc%p = reshape([B_V_c(:,in_limits(1):in_limits(2)),&
                    &B_V_s(:,in_limits(1):in_limits(2))],[2*mnmax_V*n_r_eq])
                
                ! jac_V
                in_1D_loc => in_1D(id); id = id+1
                in_1D_loc%var_name = 'jac_V'
                allocate(in_1D_loc%tot_i_min(3),in_1D_loc%tot_i_max(3))
                allocate(in_1D_loc%loc_i_min(3),in_1D_loc%loc_i_max(3))
                in_1D_loc%loc_i_min = [1,1,1]
                in_1D_loc%loc_i_max = [mnmax_V,n_r_eq,2]
                in_1D_loc%tot_i_min = in_1D_loc%loc_i_min
                in_1D_loc%tot_i_max = in_1D_loc%loc_i_max
                allocate(in_1D_loc%p(2*mnmax_V*n_r_eq))
                in_1D_loc%p = reshape([jac_V_c(:,in_limits(1):in_limits(2)),&
                    &jac_V_s(:,in_limits(1):in_limits(2))],[2*mnmax_V*n_r_eq])
#endif
            case (2)                                                            ! HELENA
                ! misc_in_H
                in_1D_loc => in_1D(id); id = id+1
                in_1D_loc%var_name = 'misc_in_H'
                allocate(in_1D_loc%tot_i_min(1),&
                    &in_1D_loc%tot_i_max(1))
                allocate(in_1D_loc%loc_i_min(1),&
                    &in_1D_loc%loc_i_max(1))
                in_1D_loc%loc_i_min = [1]
                in_1D_loc%loc_i_max = [2]
                in_1D_loc%tot_i_min = in_1D_loc%loc_i_min
                in_1D_loc%tot_i_max = in_1D_loc%loc_i_max
                allocate(in_1D_loc%p(2))
                in_1D_loc%p = [ias*1._dp,nchi*1._dp]
                
                ! RZ_H
                in_1D_loc => in_1D(id); id = id+1
                in_1D_loc%var_name = 'RZ_H'
                allocate(in_1D_loc%tot_i_min(3),in_1D_loc%tot_i_max(3))
                allocate(in_1D_loc%loc_i_min(3),in_1D_loc%loc_i_max(3))
                in_1D_loc%loc_i_min = [1,1,1]
                in_1D_loc%loc_i_max = [nchi,n_r_eq,2]
                in_1D_loc%tot_i_min = in_1D_loc%loc_i_min
                in_1D_loc%tot_i_max = in_1D_loc%loc_i_max
                allocate(in_1D_loc%p(2*nchi*n_r_eq))
                in_1D_loc%p = reshape([R_H(:,in_limits(1):in_limits(2)),&
                    &Z_H(:,in_limits(1):in_limits(2))],[2*nchi*n_r_eq])
                
                ! chi_H
                in_1D_loc => in_1D(id); id = id+1
                in_1D_loc%var_name = 'chi_H'
                allocate(in_1D_loc%tot_i_min(1),in_1D_loc%tot_i_max(1))
                allocate(in_1D_loc%loc_i_min(1),in_1D_loc%loc_i_max(1))
                in_1D_loc%loc_i_min = [1]
                in_1D_loc%loc_i_max = [nchi]
                in_1D_loc%tot_i_min = in_1D_loc%loc_i_min
                in_1D_loc%tot_i_max = in_1D_loc%loc_i_max
                allocate(in_1D_loc%p(nchi))
                in_1D_loc%p = chi_H
                
                ! flux_p_H
                in_1D_loc => in_1D(id); id = id+1
                in_1D_loc%var_name = 'flux_p_H'
                allocate(in_1D_loc%tot_i_min(1),in_1D_loc%tot_i_max(1))
                allocate(in_1D_loc%loc_i_min(1),in_1D_loc%loc_i_max(1))
                in_1D_loc%loc_i_min = [1]
                in_1D_loc%loc_i_max = [n_r_eq]
                in_1D_loc%tot_i_min = in_1D_loc%loc_i_min
                in_1D_loc%tot_i_max = in_1D_loc%loc_i_max
                allocate(in_1D_loc%p(n_r_eq))
                in_1D_loc%p = flux_p_H(in_limits(1):in_limits(2))
                
                ! qs_H
                in_1D_loc => in_1D(id); id = id+1
                in_1D_loc%var_name = 'qs_H'
                allocate(in_1D_loc%tot_i_min(1),in_1D_loc%tot_i_max(1))
                allocate(in_1D_loc%loc_i_min(1),in_1D_loc%loc_i_max(1))
                in_1D_loc%loc_i_min = [1]
                in_1D_loc%loc_i_max = [n_r_eq]
                in_1D_loc%tot_i_min = in_1D_loc%loc_i_min
                in_1D_loc%tot_i_max = in_1D_loc%loc_i_max
                allocate(in_1D_loc%p(n_r_eq))
                in_1D_loc%p = qs_H(in_limits(1):in_limits(2))
                
                ! pres_H
                in_1D_loc => in_1D(id); id = id+1
                in_1D_loc%var_name = 'pres_H'
                allocate(in_1D_loc%tot_i_min(1),in_1D_loc%tot_i_max(1))
                allocate(in_1D_loc%loc_i_min(1),in_1D_loc%loc_i_max(1))
                in_1D_loc%loc_i_min = [1]
                in_1D_loc%loc_i_max = [n_r_eq]
                in_1D_loc%tot_i_min = in_1D_loc%loc_i_min
                in_1D_loc%tot_i_max = in_1D_loc%loc_i_max
                allocate(in_1D_loc%p(n_r_eq))
                in_1D_loc%p = pres_H(in_limits(1):in_limits(2))
                
                ! RBphi_H
                in_1D_loc => in_1D(id); id = id+1
                in_1D_loc%var_name = 'RBphi_H'
                allocate(in_1D_loc%tot_i_min(1),in_1D_loc%tot_i_max(1))
                allocate(in_1D_loc%loc_i_min(1),in_1D_loc%loc_i_max(1))
                in_1D_loc%loc_i_min = [1]
                in_1D_loc%loc_i_max = [n_r_eq]
                in_1D_loc%tot_i_min = in_1D_loc%loc_i_min
                in_1D_loc%tot_i_max = in_1D_loc%loc_i_max
                allocate(in_1D_loc%p(n_r_eq))
                in_1D_loc%p = RBphi_H(in_limits(1):in_limits(2))
                
#if ldebug
                ! h_H
                in_1D_loc => in_1D(id); id = id+1
                in_1D_loc%var_name = 'h_H'
                allocate(in_1D_loc%tot_i_min(3),in_1D_loc%tot_i_max(3))
                allocate(in_1D_loc%loc_i_min(3),in_1D_loc%loc_i_max(3))
                in_1D_loc%loc_i_min = [1,1,1]
                in_1D_loc%loc_i_max = [nchi,n_r_eq,3]
                in_1D_loc%tot_i_min = in_1D_loc%loc_i_min
                in_1D_loc%tot_i_max = in_1D_loc%loc_i_max
                allocate(in_1D_loc%p(3*nchi*n_r_eq))
                in_1D_loc%p = reshape([h_H_11(:,in_limits(1):in_limits(2)),&
                    &h_H_12(:,in_limits(1):in_limits(2)),&
                    &h_H_33(:,in_limits(1):in_limits(2))],&
                    &[3*nchi*n_r_eq])
#endif
            case default
                err_msg = 'No equilibrium style associated with '//&
                    &trim(i2str(eq_style))
                ierr = 1
                CHCKERR(err_msg)
        end select
        
        ! misc_X
        in_1D_loc => in_1D(id); id = id+1
        in_1D_loc%var_name = 'misc_X'
        allocate(in_1D_loc%tot_i_min(1),in_1D_loc%tot_i_max(1))
        allocate(in_1D_loc%loc_i_min(1),in_1D_loc%loc_i_max(1))
        in_1D_loc%loc_i_min = [1]
        in_1D_loc%loc_i_max = [9]
        in_1D_loc%tot_i_min = in_1D_loc%loc_i_min
        in_1D_loc%tot_i_max = in_1D_loc%loc_i_max
        allocate(in_1D_loc%p(9))
        in_1D_loc%p = [prim_X*1._dp,n_mod_X*1._dp,min_sec_X*1._dp,&
            &max_sec_X*1._dp,norm_disc_prec_X*1._dp,norm_style*1._dp,&
            &U_style*1._dp,X_style*1._dp,matrix_SLEPC_style*1._dp]
        
        ! misc_sol
        in_1D_loc => in_1D(id); id = id+1
        in_1D_loc%var_name = 'misc_sol'
        allocate(in_1D_loc%tot_i_min(1),in_1D_loc%tot_i_max(1))
        allocate(in_1D_loc%loc_i_min(1),in_1D_loc%loc_i_max(1))
        in_1D_loc%loc_i_min = [1]
        in_1D_loc%loc_i_max = [7]
        in_1D_loc%tot_i_min = in_1D_loc%loc_i_min
        in_1D_loc%tot_i_max = in_1D_loc%loc_i_max
        allocate(in_1D_loc%p(7))
        in_1D_loc%p = [min_r_sol,max_r_sol,alpha,&
            &norm_disc_prec_sol*1._dp,BC_style*1._dp,EV_style*1._dp,EV_BC]
        
        call lvl_ud(-1)
        
        ! write
        call writo('Writing using HDF5')
        call lvl_ud(1)
        ierr = print_HDF5_arrs(in_1D(1:id-1),PB3D_name,'in')
        CHCKERR('')
        call lvl_ud(-1)
        
        ! deallocate
        deallocate(in_1D)
        
        ! clean up
        nullify(in_1D_loc)
        
        ! user output
        call lvl_ud(-1)
        call writo('Input variables written to output')
    end function print_output_in
    
    ! Cleans up input from equilibrium codes.
    integer function dealloc_in() result(ierr)
        use num_vars, only: eq_style
        use VMEC, only: dealloc_VMEC
        use HELENA_vars, only: dealloc_HEL
        
        character(*), parameter :: rout_name = 'dealloc_in'
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! deallocate depending on equilibrium style
        select case (eq_style)
            case (1)                                                            ! VMEC
                call dealloc_VMEC
            case (2)                                                            ! HELENA
                call dealloc_HEL
            case default
                err_msg = 'No equilibrium style associated with '//&
                    &trim(i2str(eq_style))
                ierr = 1
                CHCKERR(err_msg)
        end select
    end function dealloc_in
    
    ! closes the output file
    ! [MPI] only master
    subroutine close_output
        use num_vars, only: rank, output_i
        
        call writo('Closing output files')
        call writo('')
        
        if (rank.eq.0) close(output_i)
    end subroutine
end module files_ops
