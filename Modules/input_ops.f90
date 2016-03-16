!------------------------------------------------------------------------------!
!   Operations concerning giving input                                         !
!------------------------------------------------------------------------------!
module input_ops
#include <PB3D_macros.h>
    use str_ops
    use output_ops
    use messages
    use num_vars, only: dp, max_str_ln
    
    implicit none
    private
    public read_input_opts, read_input_eq, print_output_in
    
contains
    ! reads input options from user-provided input file
    ! These variables will be passed on through MPI.
    ! [MPI] only master
    integer function read_input_opts() result(ierr)
        use num_vars, only: &
            &max_it_NR, tol_NR, max_it_rich, input_i, use_pol_flux_F, &
            &EV_style, max_mem_per_proc, plot_resonance, tol_rich, &
            &n_sol_requested, rank, plot_grid, plot_flux_q, &
            &use_normalization, n_sol_plotted, n_theta_plot, n_zeta_plot, &
            &EV_BC, rho_style, retain_all_sol, prog_style, norm_disc_prec_X, &
            &norm_disc_prec_eq, norm_disc_prec_sol, BC_style, max_it_inv, &
            &tol_norm, tol_SLEPC, max_it_slepc, n_procs, pi, plot_size, &
            &U_style, norm_style, test_max_mem, X_style, matrix_SLEPC_style,  &
            &input_name, rich_restart_lvl, eq_style
        use eq_vars, only: rho_0, R_0, pres_0, B_0, psi_0, T_0
        use messages, only: writo, lvl_ud
        use X_vars, only: min_r_sol, max_r_sol, n_mod_X, prim_X, min_sec_X, &
            &max_sec_X
        use rich_vars, only: rich_lvl, min_n_par_X, req_min_n_par_X
        use rich_ops, only: find_max_lvl_rich
        use grid_vars, only: n_r_sol, min_par_X, max_par_X
        use sol_vars, only: alpha
        
        character(*), parameter :: rout_name = 'read_input_opts'
        
        ! local variables
        integer :: PB3D_lvl_rich                                                ! Richardson level to post-process (for POST)
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! input options
        namelist /inputdata_PB3D/ min_par_X, max_par_X, alpha, min_r_sol, &
            &max_r_sol, max_it_NR, tol_NR, use_pol_flux_F, rho_style, &
            &rho_0, plot_grid, plot_flux_q, prim_X, min_sec_X, max_sec_X, &
            &n_mod_X, use_normalization, n_theta_plot, n_zeta_plot, &
            &norm_disc_prec_eq, tol_norm, max_mem_per_proc, n_r_sol, &
            &max_it_rich, tol_rich, EV_style, plot_resonance, n_sol_requested, &
            &EV_BC, tol_SLEPC, retain_all_sol, pres_0, R_0, psi_0, B_0, T_0, &
            &norm_disc_prec_X, BC_style, max_it_inv, max_it_slepc, &
            &norm_disc_prec_sol, plot_size, U_style, norm_style, test_max_mem, &
            &matrix_SLEPC_style, rich_restart_lvl, min_n_par_X
        namelist /inputdata_POST/ n_sol_plotted, n_theta_plot, n_zeta_plot, &
            &plot_resonance, plot_flux_q, plot_grid, norm_disc_prec_sol, &
            &plot_size, test_max_mem, PB3D_lvl_rich, max_it_NR, tol_NR
        
        ! initialize ierr
        ierr = 0
        
        if (rank.eq.0) then                                                     ! only master
            call writo('Setting up user-provided input "'//trim(input_name)//&
                &'"')
            call lvl_ud(1)
            
            ! initialize input variables (optionally overwritten by user later)
            call writo('Initialize all the inputs with default values')
            call lvl_ud(1)
            
            ! common variables for all program styles
            max_mem_per_proc = 6000_dp/n_procs                                  ! count with 6GB
            plot_size = [10,5]                                                  ! size of plot in inch
            test_max_mem = .false.                                              ! do not test maximum memory
            select case(eq_style)
                case (1)                                                        ! VMEC
                    n_theta_plot = 201                                          ! nr. poloidal points in plot
                    n_zeta_plot = 101                                           ! nr. toroidal points in plot
                case (2)                                                        ! HELENA
                    n_theta_plot = 501                                          ! nr. poloidal points in plot
                    n_zeta_plot = 1                                             ! nr. toroidal points in plot
                case default
                    err_msg = 'No equilibrium style associated with '//&
                        &trim(i2str(eq_style))
                    ierr = 1
                    CHCKERR(err_msg)
            end select
            
            ! select depending on program style
            select case (prog_style)
                case(1)                                                         ! PB3D
                    call default_input_PB3D
                case(2)                                                         ! POST
                    ierr = default_input_POST()
                    CHCKERR('')
                case default
                    err_msg = 'No program style associated with '//&
                        &trim(i2str(prog_style))
                    ierr = 1
                    CHCKERR(err_msg)
            end select
            
            call lvl_ud(-1)
            
            ! read user input
            ! select depending on program style
            select case (prog_style)
                case(1)                                                         ! PB3D
                    read(input_i,nml=inputdata_PB3D,iostat=ierr)                ! read input data
                    
                    ! check input if successful read
                    if (ierr.eq.0) then                                         ! input file succesfully read
                        call writo('Overwriting with user-provided file "'&
                            &//trim(input_name) // '"')
                        
                        call writo('Checking user-provided file')
                        
                        call lvl_ud(1)
                        
                        ! adapt run-time variables if needed
                        call adapt_run
                        
                        ! adapt plotting variables if needed
                        call adapt_plot
                        
                        ! adapt perturbation modes
                        ierr = adapt_X_modes()
                        CHCKERR('')
                        
                        ! adapt min_n_par_X if needed
                        call adapt_min_n_par_X
                        
                        ! adapt tolerances if needed
                        call adapt_tol_rich
                        
                        ! adapt solution grid
                        ierr = adapt_sol_grid()
                        CHCKERR('')
                        
                        ! adapt Newton-Rhapson variables if needed
                        call adapt_NR
                        
                        ! adapt normalization variables if needed
                        ierr = adapt_normalization()
                        CHCKERR('')
                        
                        ! adapt input / output variables if needed
                        ierr = adapt_inoutput()
                        CHCKERR('')
                        
                        ! adapt Richardson variables if needed
                        call adapt_r
                        
                        ! adapt variables for inverse if needed
                        call adapt_inv
                        
                        ! adapt Newton-Rhapson variables if needed
                        call adapt_NR
                        
                        call lvl_ud(-1)
                    else                                                        ! cannot read input data
                        ierr = 1
                        err_msg = 'Cannot read user-provided file "'&
                            &//trim(input_name)//'"'
                        CHCKERR(err_msg)
                    end if
                    
                    ! multiply alpha by pi
                    alpha = alpha*pi
                case(2)                                                         ! POST
                    read(input_i,nml=inputdata_POST,iostat=ierr)                ! read input data
                    
                    ! check input if successful read
                    if (ierr.eq.0) then                                         ! input file succesfully read
                        call writo('Overwriting with user-provided file "'&
                            &//trim(input_name) // '"')
                    else                                                        ! cannot read input data
                        ierr = 1
                        err_msg = 'Cannot read user-provided file "'&
                            &//trim(input_name)//'"'
                        CHCKERR(err_msg)
                    end if
                    
                    ! set up max_it_rich and rich_lvl
                    max_it_rich = PB3D_lvl_rich
                    rich_lvl = PB3D_lvl_rich
                case default
                    err_msg = 'No program style associated with '//&
                        &trim(i2str(prog_style))
                    ierr = 1
                    CHCKERR(err_msg)
            end select
            
            ! user output
            call writo('Close input file')
            close(input_i)
            
            call lvl_ud(-1)
            call writo('Input values set')
        end if
    contains
        subroutine default_input_PB3D
            use num_vars, only: use_pol_flux_E
            
            ! concerning Newton-Rhapson
            max_it_NR = 500                                                     ! maximum 500 Newton-Rhapson iterations
            tol_NR = 1.0E-10_dp                                                 ! wanted relative error in Newton-Rhapson iteration
            
            ! runtime variables
            use_normalization = .true.                                          ! use normalization for the variables
            norm_disc_prec_eq = 1                                               ! precision 1 normal discretization of equilibrium
            norm_disc_prec_X = 1                                                ! precision 1 normal discretization of perturbation
            norm_disc_prec_sol = 1                                              ! precision 1 normal discretization of solution
            max_it_slepc = 10000                                                ! max. nr. of iterations for SLEPC
            EV_BC = 1._dp                                                       ! use 1 as artificial EV for the Boundary Conditions
            tol_SLEPC = 1.E-8_dp                                                ! tolerance of 1E-8
            rho_style = 1                                                       ! constant pressure profile, equal to rho_0
            U_style = 3                                                         ! full expression for U, up to order 3
            norm_style = 1                                                      ! perpendicular kinetic energy normalized
            BC_style = [1,2]                                                    ! left BC zeroed and right BC through minimization of energy
            X_style = 2                                                         ! fast style: mode numbers optimized in normal coordinate
            matrix_SLEPC_style = 1                                              ! sparse matrix storage
            plot_resonance = .false.                                            ! do not plot the q-profile with nq-m = 0
            plot_grid = .false.                                                 ! do not plot the grid
            plot_flux_q = .false.                                               ! do not plot the flux quantities
            
            ! variables concerning input / output
            n_sol_requested = 3                                                 ! request solutions with 3 highes EV
            retain_all_sol = .false.                                            ! don't retain faulty ones
            rich_restart_lvl = 1                                                ! don't restart
            
            ! variables concerning poloidal mode numbers m
            min_n_par_X = req_min_n_par_X                                       ! nonsensible value to check for user overwriting
            prim_X = 20                                                         ! main mode number of perturbation
            min_sec_X = huge(1)                                                 ! nonsensible value to check for user overwriting
            max_sec_X = huge(1)                                                 ! nonsensible value to check for user overwriting
            n_mod_X = huge(1)                                                   ! nonsensible value to check for user overwriting
            min_par_X = -4.0_dp                                                 ! minimum parallel angle [pi]
            max_par_X = 4.0_dp                                                  ! maximum parallel angle [pi]
            use_pol_flux_F = use_pol_flux_E                                     ! use same normal flux coordinate as the equilibrium
            
            ! variables concerning alpha
            alpha = 0._dp                                                       ! field line based in outboard
            
            ! variables concerning solution
            min_r_sol = 0.1_dp                                                  ! minimum normal range
            max_r_sol = 1.0_dp                                                  ! maximum normal range
            tol_norm = 0.05                                                     ! tolerance for normal range
            EV_style = 1                                                        ! slepc solver for EV problem
            n_r_sol = 100                                                       ! 100 points in solution grid
            
            ! variables concerning normalization
            rho_0 = 10E-6_dp                                                    ! for fusion, particle density of around 1E21, mp around 1E-27
            R_0 = huge(1._dp)                                                   ! nonsensible value to check for user overwriting
            pres_0 = huge(1._dp)                                                ! nonsensible value to check for user overwriting
            psi_0 = huge(1._dp)                                                 ! nonsensible value to check for user overwriting
            B_0 = huge(1._dp)                                                   ! nonsensible value to check for user overwriting
            T_0 = huge(1._dp)                                                   ! nonsensible value to check for user overwriting
            
            ! concerning Richardson extrapolation
            max_it_rich = 1                                                     ! by default no Richardson extrapolation
            tol_rich = 1E-5                                                     ! wanted relative error in Richardson extrapolation
            
            ! concerning calculating the inverse
            max_it_inv = 1                                                      ! by default no iteration to calculate inverse
        end subroutine default_input_PB3D
        
        integer function default_input_POST() result(ierr)
            character(*), parameter :: rout_name = 'default_input_POST'
            
            ! initialize ierr
            ierr = 0
            
            ! concerning Newton-Rhapson
            max_it_NR = 5000                                                    ! more iterations than PB3D
            tol_NR = 1.0E-8_dp                                                  ! less relative error than PB3D
            
            ! runtime variables
            plot_resonance = .false.                                            ! do not plot the q-profile with nq-m = 0
            plot_flux_q = .false.                                               ! do not plot the flux quantities
            plot_grid = .false.                                                 ! do not plot the grid
            norm_disc_prec_sol = 1                                              ! precision 1 normal discretization of solution
            
            ! Richardson variables
            ierr = find_max_lvl_rich(PB3D_lvl_rich)                             ! highest Richardson level found
            CHCKERR('')
            if (PB3D_lvl_rich.lt.0) then
                ierr = 1
                err_msg = 'No suitable Richardson level found'
                CHCKERR(err_msg)
            end if
            call writo('Maximum Richardson level found: '//&
                &trim(i2str(PB3D_lvl_rich)))
            
            ! variables concerning output
            n_sol_plotted = n_sol_requested                                     ! plot all solutions
        end function default_input_POST
        
        ! checks whether the variables concerning run-time are chosen correctly.
        ! rho_style has  to be 1  (constant rho = rho_0)  and matrix_SLEPC_style
        ! has to be 0..1.
        subroutine adapt_run
            use num_vars, only: eq_style
            
            if (rho_style.ne.1) then
                ierr = 1
                err_msg = 'rho_style has to be 1 (constant)'
                CHCKERR(err_msg)
            end if
            if (eq_style.eq.2 .and. .not.use_normalization) then
                use_normalization = .true.
                call writo('WARNING: normalization is always used for HELENA')
            end if
            if (matrix_SLEPC_style.lt.1 .or. matrix_SLEPC_style.gt.2) then
                ierr = 1
                err_msg = 'matrix_SLEPC_style has to be 1 (sparse) or 2 (shell)'
                CHCKERR(err_msg)
            end if
        end subroutine adapt_run
        
        ! checks whether the variables concerning output are chosen correctly.
        ! n_sol_requested has to be at least one
        ! rich_restart_lvl has to not more than max_it_rich
        integer function adapt_inoutput() result(ierr)
            character(*), parameter :: rout_name = 'adapt_inoutput'
            
            ! initialize ierr
            ierr = 0
            
            if (n_sol_requested.lt.1) then
                n_sol_requested = 1
                call writo('WARNING: n_sol_requested has been increased to '&
                    &//trim(i2str(n_sol_requested)))
            end if
            if (rich_restart_lvl.lt.1) then
                call writo('WARNING: rich_restart_lvl was '//&
                    &trim(i2str(rich_restart_lvl))//&
                    &' but cannot be lower than 1, so it was reset to 1')
                rich_restart_lvl = 1
            end if
            if (rich_restart_lvl.gt.max_it_rich) then
                ierr = 1
                err_msg = 'rich_restart_lvl not within 1..max_it_rich = 1..'//&
                    &trim(i2str(max_it_rich))
                CHCKERR(err_msg)
            end if
        end function adapt_inoutput
        
        ! checks whether the variables concerning plotting are chosen correctly.
        ! n_theta and n_zeta_plot have to be positive
        subroutine adapt_plot
            if (n_theta_plot.lt.1) then
                n_theta_plot = 1
                call writo('WARNING: n_theta_plot cannot be negative and &
                    &is set to '//trim(i2str(n_theta_plot)))
            end if
            if (n_zeta_plot.lt.1) then
                n_zeta_plot = 1
                call writo('WARNING: n_zeta_plot cannot be negative and &
                    &is set to '//trim(i2str(n_zeta_plot)))
            end if
        end subroutine adapt_plot
        
        ! checks whether  variables concerning  perturbation modes  are correct:
        ! the absolute value of prim_X has to be at least min_nm_X.
        ! Also  deduces the X  style; min_sec_X  and max_sec_X, pertaining  to X
        ! style 1 and n_sec_X, pertaining to X style 2 cannot be both given.
        ! Furthermore, if  X style 1  (prescribed), max_sec_X has to  be greater
        ! than min_sec_X in absolute value. If  X style 2 (fast), n_mod_X has to
        ! be a number greater than 0.
        ! Lastly,  for every X  style, some checks  are made and  some variables
        ! set.
        integer function adapt_X_modes() result(ierr)
            use X_vars, only: min_nm_X
            
            character(*), parameter :: rout_name = 'adapt_X_modes'
            
            ! local variables
            character(len=max_str_ln) :: err_msg                                ! error message
            
            ! initialize ierr
            ierr = 0
            
            ! check prim_X: common for all X styles
            if (abs(prim_X).lt.min_nm_X) then
                ierr = 1
                err_msg = 'prim_X should be at least '//trim(i2str(min_nm_X))//&
                    &', and probably a bit higher still'
                CHCKERR(err_msg)
            end if
            
            ! set X_style
            if (min_sec_X.ne.huge(1) .and. min_sec_X.ne.huge(1)) then           ! min_sec_X and max_sec_X prescribed
                if (n_mod_X.eq.huge(1)) then                                    ! n_mod_X not prescribed
                    X_style = 1                                                 ! prescribed style
                else                                                            ! n_mod_X also prescribed
                    ierr = 1
                    err_msg = 'Cannot prescribe both min and max of secondary &
                        &X and nr. of modes.'
                    CHCKERR(err_msg)
                end if
            else                                                                ! min_sec_X and max_sec_X not prescribed
                if (n_mod_X.eq.huge(1)) then                                    ! n_mod_X not prescribed either: default
                    X_style = 2                                                 ! X style 2: fast
                    n_mod_X = 10                                                ! 10 resonating modes at every place
                else                                                            ! but n_mod_X prescribed
                    X_style = 2
                end if
            end if
            
            ! checks for each X style
            select case (X_style)
                case (1)                                                        ! prescribed
                    ! check min_sec_X, max_sec_X
                    if (abs(max_sec_X).lt.abs(min_sec_X)) then
                        ierr = 1
                        err_msg = 'max_sec_X has to be larger or equal to &
                            &min_sec_X in absolute value'
                        CHCKERR(err_msg)
                    end if
                    
                    ! set n_mod_X
                    n_mod_X = abs(max_sec_X)-(min_sec_X)+1
                case (2)                                                        ! fast
                    if (n_mod_X.lt.1) then
                        ierr = 1
                        err_msg = 'the number of resonating modes has to be at &
                            &least one'
                        CHCKERR(err_msg)
                    end if
                case default
                    err_msg = 'No X style associated with '//&
                        &trim(i2str(X_style))
                    ierr = 1
                    CHCKERR(err_msg)
            end select
        end function adapt_X_modes
        
        ! Checks whether  min_n_par_X has  been chosen  correctly: It  should be
        ! larger than 10 (arbitrary) if an integral is to be calculated.
        subroutine adapt_min_n_par_X
            if (min_n_par_X.lt.req_min_n_par_X)  then
                min_n_par_X = req_min_n_par_X
                call writo('WARNING: min_n_par_X has been increased to '//&
                    &trim(i2str(min_n_par_X)))
            end if
        end subroutine adapt_min_n_par_X
        
        ! Checks  whether variables  concerning the  solution grid  are correct:
        ! min_r_sol  should not  be too  close to  zero because  the equilibrium
        ! calculations yield an infinity at  the magnetic axis. max_r_sol cannot
        ! be larger  than 1.0  and has  to be larger  than min_r_sol.  Also, the
        ! number of normal points has to be big enough.
        integer function adapt_sol_grid() result(ierr)
            use grid_vars, only: n_r_eq
            
            character(*), parameter :: rout_name = 'adapt_sol_grid'
            
            ! local variables
            character(len=max_str_ln) :: err_msg                                ! error message
            real(dp) :: one = 1.00000001_dp                                     ! one plus a little margin
            
            ! initialize ierr
            ierr = 0
            
            ! check min_r_sol
            if (min_r_sol.lt.one/(n_r_eq-1)) then
                min_r_sol = one/(n_r_eq-1)
                call writo('WARNING: min_r_sol has been increased to '//&
                    &trim(r2strt(min_r_sol)))
            end if
            
            ! check if max_r_sol is not greater than 1
            if (max_r_sol.gt.1.0) then
                max_r_sol = 1.
                call writo('WARNING: max_r_sol has been decreased to '//&
                    &trim(r2strt(max_r_sol)))
            end if
            
            ! check if min_r_sol < max_r_sol  with at least one equilbrium point
            ! between them
            if (min_r_sol+1./(n_r_eq-1).ge.max_r_sol) then
                ierr = 1
                err_msg = 'max_r_sol - min_r_sol has to be at least '//&
                    &trim(r2strt(1._dp/(n_r_eq-1)))
                CHCKERR(err_msg)
            end if
            
            ! check n_r_sol
            if (n_r_sol.lt.6*norm_disc_prec_X+2) then
                n_r_sol = 6*norm_disc_prec_X+2
                call writo('WARNING: n_r_sol has been increased to '//&
                    &trim(i2str(n_r_sol)))
            end if
        end function adapt_sol_grid
        
        ! checks  whether the variables concerning  Richardson extrapolation are
        ! correct. max_it_rich has to be at least 1
        subroutine adapt_r
            if (max_it_rich.lt.1) then
                max_it_rich = 1
                call writo('WARNING: max_it_rich has been increased to 1')
            end if
        end subroutine adapt_r
        
        ! checks whether  the variables  concerning calculating the  inverse are
        ! correct. max_it_inv has to be at least 1
        subroutine adapt_inv
            if (max_it_inv.lt.1) then
                max_it_inv = 1
                call writo('WARNING: max_it_inv has been increased to 1')
            end if
        end subroutine adapt_inv
        
        ! checks  whether the  variables concerning Newton-Rhapson  are correct.
        ! max_it_NR has to be at least 2
        subroutine adapt_NR
            if (max_it_NR.lt.1) then
                max_it_NR = 2
                call writo('WARNING: max_it_NR has been increased to 2')
            end if
        end subroutine adapt_NR
        
        ! checks whether Richardson tolerances are correct
        subroutine adapt_tol_rich
            if (tol_norm.lt.0) then
                call writo('WARNING: tol_norm has been increased to 0')
                tol_norm = 0._dp
            end if
            if (tol_norm.gt.1) then
                call writo('WARNING: tol_norm has been decreased to 1')
                tol_norm = 1._dp
            end if
        end subroutine adapt_tol_rich
        
        ! checks whether normalization variables are chosen correctly. rho_0 has
        ! to be positive
        integer function adapt_normalization() result(ierr)
            character(*), parameter :: rout_name = 'adapt_normalization'
            
            ! local variables
            character(len=max_str_ln) :: err_msg                                ! error message
            
            ! initialize ierr
            ierr = 0
            
            if (rho_0.le.0) then
                ierr = 1
                err_msg = 'rho_0 has to be positive'
                CHCKERR(err_msg)
            end if
        end function adapt_normalization
    end function read_input_opts
    
    ! Reads the equilibrium input file if no Richardson restart.
    ! These variables will be passed on through HDF5 data files.
    integer function read_input_eq() result(ierr)
        use num_vars, only: eq_style, use_pol_flux_E, eq_i
        use VMEC, only: read_VMEC
        use HELENA_vars, only: read_HEL
        use grid_vars, only: n_r_in
        
        character(*), parameter :: rout_name = 'read_input_eq'
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! choose which equilibrium style is being used:
        !   1:  VMEC
        !   2:  HELENA
        select case (eq_style)
            case (1)                                                            ! VMEC
                ierr = read_VMEC(n_r_in,use_pol_flux_E)
                CHCKERR('')
            case (2)                                                            ! HELENA
                ierr = read_HEL(n_r_in,use_pol_flux_E)
                CHCKERR('')
            case default
                err_msg = 'No equilibrium style associated with '//&
                    &trim(i2str(eq_style))
                ierr = 1
                CHCKERR(err_msg)
        end select
        
        ! close equilibrium file
        close(eq_i)
    end function read_input_eq
    
    ! Print input quantities to an output file:
    !   - misc_in:   prog_version, eq_style, rho_style, R_0, pres_0, B_0, psi_0
    !                rho_0, T_0, vac_perm, use_pol_flux_F, use_pol_flux_E,
    !                use_normalization, norm_disc_prec_eq, n_r_in, n_r_eq,
    !                n_r_sol
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
        use grid_vars, onLy: n_r_in, n_r_eq, n_r_sol
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
        
        ! calculate limits of input range
        ierr = calc_norm_range(in_limits=in_limits)
        CHCKERR('')
        n_r_eq = in_limits(2)-in_limits(1)+1
        call writo('Only the normal range '//trim(i2str(in_limits(1)))//'..'//&
            &trim(i2str(in_limits(2)))//' is written')
        call writo('(relevant to solution range '//&
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
            &n_r_sol*1._dp,max_flux_E,max_flux_F]
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
        
        ! write
        ierr = print_HDF5_arrs(in_1D(1:id-1),PB3D_name,'in')
        CHCKERR('')
        
        ! deallocate
        deallocate(in_1D)
        
        ! clean up
        nullify(in_1D_loc)
        
        ! user output
        call lvl_ud(-1)
        call writo('Input variables written to output')
    end function print_output_in
end module input_ops
