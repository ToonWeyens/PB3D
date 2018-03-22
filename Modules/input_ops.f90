!------------------------------------------------------------------------------!
!> Operations concerning giving input.
!------------------------------------------------------------------------------!
module input_ops
#include <PB3D_macros.h>
    use str_utilities
    use output_ops
    use messages
    use num_vars, only: dp, max_str_ln
    
    implicit none
    private
    public read_input_opts, read_input_eq, print_output_in
    
contains
    !> Reads input options from user-provided input file.
    !! 
    !! This  is done  by the  master process  and these  variables will  then be
    !! passed on through MPI to other processes in broadcast_input_opts().
    !!
    !! \return ierr
    integer function read_input_opts() result(ierr)
        use num_vars, only: &
            &max_it_zero, tol_zero, max_it_rich, input_i, use_pol_flux_F, &
            &EV_style, max_tot_mem, plot_resonance, tol_rich, &
            &n_sol_requested, rank, plot_magn_grid, plot_B, plot_J, &
            &plot_flux_q, plot_kappa, plot_sol_xi, plot_sol_Q, plot_E_rec, &
            &use_normalization, n_sol_plotted, n_theta_plot, n_zeta_plot, &
            &EV_BC, rho_style, retain_all_sol, prog_style, norm_disc_prec_X, &
            &norm_disc_prec_eq, norm_disc_prec_sol, norm_disc_style_sol, &
            &BC_style, tol_norm, tol_SLEPC_loc => tol_SLEPC, max_it_SLEPC, &
            &plot_size, U_style, norm_style, X_style, matrix_SLEPC_style, &
            &input_name, rich_restart_lvl, eq_style, relax_fac_HH, &
            &min_theta_plot, max_theta_plot, min_zeta_plot, max_zeta_plot, &
            &min_r_plot, max_r_plot, max_nr_backtracks_HH, POST_style, &
            &plot_grid_style, def_relax_fac_HH, magn_int_style, K_style, &
            &ex_plot_style, pert_mult_factor_POST, sol_n_procs, n_procs, &
            &POST_output_full, POST_output_sol, EV_guess, solver_SLEPC_style, &
            &plot_vac_pot, min_Rvac_plot, max_Rvac_plot, min_Zvac_plot, &
            &max_Zvac_plot, n_vac_plot, alpha_style, X_grid_style, V_interp_style
        use eq_vars, only: rho_0, R_0, pres_0, B_0, psi_0, T_0
        use X_vars, only: min_r_sol, max_r_sol, n_mod_X, prim_X, min_sec_X, &
            &max_sec_X
        use rich_vars, only: rich_lvl, min_n_par_X, req_min_n_par_X
        use rich_ops, only: find_max_rich_lvl
        use grid_vars, only: n_r_sol, min_par_X, max_par_X, n_alpha, &
            &min_alpha, max_alpha
        
        character(*), parameter :: rout_name = 'read_input_opts'
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    !< error message
        integer :: PB3D_rich_lvl                                                !< Richardson level to post-process (for POST)
        integer, parameter :: max_size_tol_SLEPC = 100                          !< maximum size of tol_SLEPC
        real(dp) :: alpha                                                       !< field line label (for alpha_style 1)
        real(dp) :: tol_SLEPC(max_size_tol_SLEPC)                               !< tol_SLEPC
        real(dp), parameter :: min_tol = 1.E-12_dp                              !< minimum general tolerance
        real(dp), parameter :: max_tol = 1.E-3_dp                               !< maximum general tolerance
        
        ! input options
        namelist /inputdata_PB3D/ min_par_X, max_par_X, alpha, min_r_sol, &
            &max_r_sol, max_it_zero, tol_zero, use_pol_flux_F, rho_style, &
            &rho_0, plot_magn_grid, plot_B, plot_J, plot_flux_q, plot_kappa, &
            &prim_X, min_sec_X, max_sec_X, n_mod_X, use_normalization, &
            &n_theta_plot, n_zeta_plot, norm_disc_prec_eq, tol_norm, &
            &max_tot_mem, n_r_sol, max_it_rich, tol_rich, EV_style, EV_guess, &
            &plot_resonance, n_sol_requested, EV_BC, tol_SLEPC, sol_n_procs, &
            &retain_all_sol, pres_0, R_0, psi_0, B_0, T_0, norm_disc_prec_X, &
            &BC_style, max_it_slepc, norm_disc_prec_sol, norm_disc_style_sol, &
            &plot_size, U_style, norm_style, K_style, matrix_SLEPC_style, &
            &rich_restart_lvl, min_n_par_X, relax_fac_HH, min_theta_plot, &
            &max_theta_plot, min_zeta_plot, max_zeta_plot, alpha_style, &
            &max_nr_backtracks_HH, magn_int_style, ex_plot_style, n_alpha, &
            &solver_SLEPC_style, X_grid_style, V_interp_style
        namelist /inputdata_POST/ n_sol_plotted, n_theta_plot, n_zeta_plot, &
            &plot_resonance, plot_flux_q, plot_kappa, plot_magn_grid, plot_B, &
            &plot_J, plot_sol_xi, plot_sol_Q, plot_E_rec, norm_disc_prec_sol, &
            &plot_size, PB3D_rich_lvl, max_it_zero, tol_zero, relax_fac_HH, &
            &min_theta_plot, max_theta_plot, min_zeta_plot, max_zeta_plot, &
            &min_r_plot, max_r_plot, max_nr_backtracks_HH, POST_style, &
            &plot_grid_style, max_tot_mem, ex_plot_style, plot_vac_pot, &
            &pert_mult_factor_POST, min_Rvac_plot, max_Rvac_plot, &
            &min_Zvac_plot, max_Zvac_plot, n_vac_plot
        
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
            sol_n_procs = n_procs                                               ! use equal amount of processes for solution as for rest of program
            max_tot_mem = 6000_dp                                               ! count with 6GB total
            plot_size = [10,5]                                                  ! size of plot in inch
            min_r_plot = 0._dp                                                  ! minimal plot normal range
            max_r_plot = 1._dp                                                  ! maximal plot normal range
            min_theta_plot = 1                                                  ! starting from pi gives nicer plots
            max_theta_plot = 3
            select case(eq_style)
                case (1)                                                        ! VMEC
                    n_theta_plot = 201                                          ! nr. poloidal points in plot
                    n_zeta_plot = 51                                            ! nr. toroidal points in plot
                    min_zeta_plot = 0                                           ! min. toroidal plot angle [pi]
                    max_zeta_plot = 2                                           ! max. toroidal plot angle [pi]
                case (2)                                                        ! HELENA
                    n_theta_plot = 501                                          ! nr. poloidal points in plot
                    n_zeta_plot = 1                                             ! nr. toroidal points in plot
                    min_zeta_plot = 0                                           ! min. toroidal plot angle [pi]
                    max_zeta_plot = min_zeta_plot                               ! max. toroidal plot angle [pi]
            end select
            plot_grid_style = 0                                                 ! normal plots on 3D geometry
            relax_fac_HH = def_relax_fac_HH                                     ! default relaxation factor
            max_nr_backtracks_HH = 20                                           ! standard nr. of backtracks
            ex_plot_style = 1                                                   ! GNUPlot
            
            ! select depending on program style
            select case (prog_style)
                case(1)                                                         ! PB3D
                    call default_input_PB3D()
                case(2)                                                         ! POST
                    call default_input_POST()
            end select
            
            call lvl_ud(-1)
            
            ! read user input
            ! select depending on program style
            select case (prog_style)
                case(1)                                                         ! PB3D
                    read(UNIT=input_i,NML=inputdata_PB3D,iostat=ierr,&
                        &iomsg=err_msg)                                         ! read input data
                    
                    ! check input if successful read
                    if (ierr.eq.0) then                                         ! input file succesfully read
                        call writo('Overwriting with user-provided file "'&
                            &//trim(input_name) // '"')
                        
                        call writo('Checking user-provided file')
                        
                        call lvl_ud(1)
                        
                        ! adapt MPI variables if needed
                        call adapt_MPI()
                        
                        ! adapt run-time variables if needed
                        ierr = adapt_run_PB3D()
                        CHCKERR('')
                        
                        ! adapt plotting variables if needed
                        ierr = adapt_plot()
                        CHCKERR('')
                        
                        ! adapt perturbation modes
                        ierr = adapt_X_modes()
                        CHCKERR('')
                        
                        ! adapt min_n_par_X if needed
                        call adapt_min_n_par_X
                        
                        ! adapt alpha if needed
                        ierr = adapt_alpha()
                        CHCKERR('')
                        
                        ! adapt solution grid
                        ierr = adapt_sol_grid(min_r_sol,max_r_sol,'sol')
                        CHCKERR('')
                        
                        ! adapt normalization variables if needed
                        ierr = adapt_normalization()
                        CHCKERR('')
                        
                        ! adapt Richardson variables if needed
                        call adapt_rich
                        
                        ! adapt input / output variables if needed
                        ierr = adapt_inoutput_PB3D()
                        CHCKERR('')
                        
                        ! adapt Householder variables if needed
                        call adapt_zero
                        
                        ! adapt tolerances if needed
                        call adapt_tols
                        
                        call lvl_ud(-1)
                    else                                                        ! cannot read input data
                        ierr = 1
                        call writo('Cannot read user-provided file "'&
                            &//trim(input_name)//'", error message:',&
                            &error=.true.)
                        call lvl_ud(1)
                        call writo('"'//trim(err_msg)//'"')
                        call lvl_ud(-1)
                        CHCKERR("")
                    end if
                case(2)                                                         ! POST
                    read(UNIT=input_i,NML=inputdata_POST,iostat=ierr,&
                        &iomsg=err_msg)                                         ! read input data
                    
                    ! check input if successful read
                    if (ierr.eq.0) then                                         ! input file succesfully read
                        call writo('Overwriting with user-provided file "'&
                            &//trim(input_name) // '"')
                        
                        call writo('Checking user-provided file')
                        
                        call lvl_ud(1)
                        
                        ! adapt run-time variables if needed
                        ierr = adapt_run_POST()
                        CHCKERR('')
                        
                        ! adapt plotting variables if needed
                        ierr = adapt_plot()
                        CHCKERR('')
                        
                        ! adapt plot grid and save to solution variable
                        ierr = adapt_sol_grid(min_r_plot,max_r_plot,'plot')
                        CHCKERR('')
                        
                        ! adapt input / output variables if needed
                        ierr = adapt_inoutput_POST()
                        CHCKERR('')
                        
                        call lvl_ud(-1)
                    else                                                        ! cannot read input data
                        ierr = 1
                        call writo('Cannot read user-provided file "'&
                            &//trim(input_name)//'", error message:',&
                            &error=.true.)
                        call lvl_ud(1)
                        call writo('"'//trim(err_msg)//'"')
                        call lvl_ud(-1)
                        CHCKERR("")
                    end if
            end select
            
            ! user output
            call writo('Close input file')
            close(input_i)
            
            call lvl_ud(-1)
            call writo('Input values set')
        end if
    contains
        ! PB3D version
        !> \private
        subroutine default_input_PB3D
            !use num_vars, only: use_pol_flux_E
            
            ! concerning solution
            n_r_sol = 100                                                       ! 100 points in solution grid
            min_r_sol = 0.1_dp                                                  ! minimum normal range
            max_r_sol = 1.0_dp                                                  ! maximum normal range
            tol_norm = 0.05                                                     ! tolerance for normal range
            EV_style = 1                                                        ! slepc solver for EV problem
            
            ! concerning field line
            alpha = 0._dp                                                       ! field line based in outboard
            min_alpha = 0.0_dp                                                  ! minimum field-line label [pi]
            max_alpha = 2.0_dp                                                  ! maximum field-line label [pi]
            select case (eq_style)
                case (1)                                                        ! VMEC
                    alpha_style = 2                                             ! multiple field-lines, single turns
                    n_alpha = 10
                case (2)                                                        ! HELENA
                    alpha_style = 1                                             ! single field-line, multiple turns
            end select
            
            ! concerning perturbation
            min_n_par_X = req_min_n_par_X                                       ! nonsensible value to check for user overwriting
            min_par_X = 0.0_dp                                                  ! minimum parallel angle [pi]
            max_par_X = 2.0_dp                                                  ! maximum parallel angle [pi]
            prim_X = 20                                                         ! main mode number of perturbation
            min_sec_X = huge(1)                                                 ! nonsensible value to check for user overwriting
            max_sec_X = huge(1)                                                 ! nonsensible value to check for user overwriting
            n_mod_X = huge(1)                                                   ! nonsensible value to check for user overwriting
            !use_pol_flux_F = use_pol_flux_E                                     ! use same normal flux coordinate as the equilibrium
            use_pol_flux_F = .true.                                             ! use poloidail flux as toroidal flux is untested
            
            ! concerning normalization
            rho_0 = 1E-7_dp                                                     ! for fusion, particle density of around 5E19 for detuerium
            R_0 = huge(1._dp)                                                   ! nonsensible value to check for user overwriting
            pres_0 = huge(1._dp)                                                ! nonsensible value to check for user overwriting
            psi_0 = huge(1._dp)                                                 ! nonsensible value to check for user overwriting
            B_0 = huge(1._dp)                                                   ! nonsensible value to check for user overwriting
            T_0 = huge(1._dp)                                                   ! nonsensible value to check for user overwriting
            
            ! concerning input / output
            n_sol_requested = 1                                                 ! request solution closest to guess
            retain_all_sol = .false.                                            ! don't retain faulty ones
            plot_resonance = .false.                                            ! do not plot the q-profile with nq-m = 0
            plot_magn_grid = .false.                                            ! do not plot the magnetic grid
            plot_kappa = .false.                                                ! do not plot the curvature
            plot_B = .false.                                                    ! do not plot the magnetic field
            plot_J = .false.                                                    ! do not plot the current
            plot_flux_q = .false.                                               ! do not plot the flux quantities
            
            ! concerning Richardson extrapolation
            max_it_rich = 1                                                     ! by default no Richardson extrapolation
            tol_rich = 1.E-4_dp                                                 ! tolerance of 1E-4
            rich_restart_lvl = 1                                                ! don't restart
            
            ! concerning finding zeros
            max_it_zero = 100                                                   ! maximum 100 iterations
            tol_zero = 1.0E-10_dp                                               ! very low tolerance for calculation of field lines
            
            ! runtime variables
            use_normalization = .true.                                          ! use normalization for the variables
            norm_disc_prec_eq = 2                                               ! precision 1 normal discretization of equilibrium
            norm_disc_prec_X = 2                                                ! precision 1 normal discretization of perturbation
            norm_disc_prec_sol = 2                                              ! precision 1 normal discretization of solution
            norm_disc_style_sol = 2                                             ! left finite differences
            magn_int_style = 1                                                  ! trapezoidal rule
            max_it_SLEPC = 1000                                                 ! max. nr. of iterations for SLEPC
            EV_BC = 1._dp                                                       ! use 1 as artificial EV for the Boundary Conditions
            EV_guess = -0.3_dp                                                  ! sensible guess
            tol_SLEPC = huge(1._dp)                                             ! nonsensible value to check for user overwriting
            rho_style = 1                                                       ! constant pressure profile, equal to rho_0
            U_style = 3                                                         ! full expression for U, up to order 3
            K_style = 1                                                         ! perpendicular kinetic energy normalized
            norm_style = 1                                                      ! MISHKA normalization
            BC_style = [1,2]                                                    ! left BC zeroed and right BC through minimization of energy
            X_style = 2                                                         ! fast style: mode numbers optimized in normal coordinate
            solver_SLEPC_style = 1                                              ! Krylov-Schur
            matrix_SLEPC_style = 1                                              ! sparse matrix storage
            X_grid_style = 1                                                    ! use equilibrium normal grid for perturbation grid as well
            V_interp_style = 1                                                  ! use finite differences
        end subroutine default_input_PB3D
        
        ! POST version
        !> \private
        subroutine default_input_POST()
            ! concerning finding zeros
            max_it_zero = 200                                                   ! more iterations than PB3D
            tol_zero = 1.0E-8_dp                                                ! less relative error than PB3D
            
            ! runtime variables
            plot_resonance = .true.                                             ! plot the q-profile with nq-m = 0
            plot_flux_q = .true.                                                ! plot the flux quantities
            plot_kappa = .true.                                                 ! plot the curvature
            plot_magn_grid = .true.                                             ! plot the magnetic grid
            plot_B = .true.                                                     ! plot the magnetic field
            plot_J = .true.                                                     ! plot the current
            plot_sol_xi = .true.                                                ! plot plasma perturbation of solution
            plot_sol_Q = .true.                                                 ! plot magnetic perturbation of solution
            plot_E_rec = .true.                                                 ! plot energy reconstruction
            plot_vac_pot = .true.                                               ! plot vacuum potential
            norm_disc_prec_sol = 1                                              ! precision 1 normal discretization of solution
            POST_style = 1                                                      ! process on extended plot grid
            write(*,*) '!!!!!! min and max need to take into account normalization'
            min_Rvac_plot = 0.1*R_0                                             ! minimum R of vacuum plot
            max_Rvac_plot = 2*R_0                                               ! maximum R of vacuum plot
            min_Zvac_plot = -R_0                                                ! minimum R of vacuum plot
            max_Zvac_plot = R_0                                                 ! maximum R of vacuum plot
            n_vac_plot = [100,100]                                              ! number of points in R and Z of vacuum plot
            
            ! variables concerning input / output
            pert_mult_factor_POST = 0._dp                                       ! factor by which to XYZ is perturbed in POST
            PB3D_rich_lvl = huge(1)                                             ! don't restart
            
            ! variables concerning output
            n_sol_plotted = -1                                                  ! plot all solutions
        end subroutine default_input_POST
        
        ! Checks  whether the  variables  concerning MPI  are chosen  correctly:
        !   sol_n_procs cannot be greater than n_procs. If it is lower than 1,
        !   the maximum amount of n_procs is used.
        !> \private
        subroutine adapt_MPI()
            use num_vars, only: n_procs
            
            if (sol_n_procs.lt.1) then
                call writo('Using the maximum number of MPI processes '//&
                    &trim(i2str(n_procs))//' for SLEPC')
                sol_n_procs = n_procs
            end if
            if (sol_n_procs.gt.n_procs) then
                call writo('Cannot use more than '//trim(i2str(n_procs))//&
                    &' MPI processes for SLEPC',warning=.true.)
                sol_n_procs = n_procs
            end if
        end subroutine adapt_MPI
        
        ! Checks whether the variables concerning run-time are chosen correctly:
        !   rho_style has  to be 1  (constant rho = rho_0),
        !   matrix_SLEPC_style has to be 0 or 1,
        !   max_it_SLEPC has to be at least 1,
        !   magnetic integral style has to be 1..2,
        !   perturbation grid style has to be 1..2,
        !   for HELENA (eq_style 1), only poloidal flux can be used.
        !> \private
        integer function adapt_run_PB3D() result(ierr)
            use num_vars, only: eq_style
            
            character(*), parameter :: rout_name = 'adapt_run_PB3D'
            
            ! initialize ierr
            ierr = 0
            
            if (rho_style.ne.1) then
                ierr = 1
                err_msg = 'rho_style has to be 1 (constant)'
                CHCKERR(err_msg)
            end if
            if (eq_style.eq.2 .and. .not.use_normalization) then
                use_normalization = .true.
                call writo('normalization is always used for HELENA',&
                    &warning=.true.)
            end if
            if (matrix_SLEPC_style.lt.1 .or. matrix_SLEPC_style.gt.2) then
                ierr = 1
                err_msg = 'matrix_SLEPC_style has to be 1 (sparse) or 2 (shell)'
                CHCKERR(err_msg)
            end if
            if (max_it_SLEPC.le.0) then
                max_it_SLEPC = 1
                call writo('max_it_SLEPC has to be positive and is set to 1',&
                    &warning=.true.)
            end if
            if (magn_int_style.lt.1 .or. magn_int_style.gt.2) then
                ierr = 1
                err_msg = 'magn_int_style has to be 1 (trapezoidal) or 2 &
                    &(Simpson 3/8 rule)'
                CHCKERR(err_msg)
            end if
            if (X_grid_style.lt.1 .or. X_grid_style.gt.2) then
                ierr = 1
                err_msg = 'X_grid_style has to be 1 (equilibrium) or 2 &
                    &(solution)'
                CHCKERR(err_msg)
            end if
            if (eq_style.eq.2 .and. .not.use_pol_flux_F) then
                use_pol_flux_F = .true.
                call writo('can only use poloidal flux for HELENA',&
                    &warning=.true.)
            end if
        end function adapt_run_PB3D
        
        ! Checks whether the variables concerning run-time are chosen correctly:
        !    POST_style  should be  1  (plotting  on extended  plot  grid) or  2
        !    (plotting on field-aligned grid also used in PB3D).
        !> \private
        integer function adapt_run_POST() result(ierr)
            character(*), parameter :: rout_name = 'adapt_run_POST'
            
            ! initialize ierr
            ierr = 0
            
            if (POST_style.lt.1 .or. POST_style.gt.2) then
                ierr = 1
                err_msg = 'POST_style has to be 1 (extended grid) or 2 &
                    &(field-aligned grid)'
                CHCKERR(err_msg)
            end if
        end function adapt_run_POST
        
        ! Checks whether the variables concerning output are chosen correctly:
        !   PB3D_rich_lvl can be at most the maximally found level
        !   POST_output_full is true if non-flux quantities are plot
        !   POST_output_sol is true if solution quantities are plot
        !> \private
        integer function adapt_inoutput_POST() result(ierr)
            use num_vars, only: compare_tor_pos
            
            character(*), parameter :: rout_name = 'adapt_inoutput_POST'
            
            ! local variables
            integer :: max_PB3D_rich_lvl                                        ! maximum Richardson level
            
            ! initialize ierr
            ierr = 0
            
            ! find maximum levl
            ierr = find_max_rich_lvl(PB3D_rich_lvl,max_PB3D_rich_lvl)
            CHCKERR('')
            if (max_PB3D_rich_lvl.le.0) then
                call writo('No suitable Richardson level found, so only &
                    &equilibrium output will be done.',alert=.true.)
                plot_sol_xi = .false.
                plot_sol_Q = .false.
                plot_E_rec = .false.
                PB3D_rich_lvl = 1
            else
                if (PB3D_rich_lvl.eq.huge(1)) then
                    PB3D_rich_lvl = max_PB3D_rich_lvl
                    call writo('Treating maximum Richardson level found: '//&
                        &trim(i2str(max_PB3D_rich_lvl)))
                else if (PB3D_rich_lvl.le.max_PB3D_rich_lvl) then
                    call writo('Treating Richardson level requested: '//&
                        &trim(i2str(PB3D_rich_lvl)))
                    if (PB3D_rich_lvl.lt.max_PB3D_rich_lvl) then
                        call lvl_ud(1)
                        call writo('Less than maximum Richardson level: '//&
                            &trim(i2str(max_PB3D_rich_lvl)),warning=.true.)
                        call lvl_ud(-1)
                    end if
                else if (PB3D_rich_lvl.gt.max_PB3D_rich_lvl) then
                    ierr = 1
                    err_msg = 'PB3D_rich_lvl cannot be higher than '//&
                        &trim(i2str(max_PB3D_rich_lvl))
                    CHCKERR(err_msg)
                else if (PB3D_rich_lvl.lt.1) then
                    call writo('Richardson level requested '//&
                        &trim(i2str(PB3D_rich_lvl))//', so only equilibrium &
                        &output will be done.',alert=.true.)
                    plot_sol_xi = .false.
                    plot_sol_Q = .false.
                    plot_E_rec = .false.
                    PB3D_rich_lvl = 1
                end if
            end if
            max_it_rich = PB3D_rich_lvl
            rich_lvl = PB3D_rich_lvl
            
            ! set POST_output_full and POST_output_sol
            POST_output_sol = plot_sol_xi .or. plot_sol_Q .or. plot_E_rec
            POST_output_full = POST_output_sol .or. plot_B .or. plot_J .or. &
                &plot_kappa .or. plot_vac_pot .or. compare_tor_pos
        end function adapt_inoutput_POST
        
        ! Checks whether the variables concerning output are chosen correctly:
        !   n_sol_requested has to be at least one,
        !   rich_restart_lvl can  be at most  one more than the  maximally found
        !   level, nor can it be higher than max_it_rich.
        !> \private
        integer function adapt_inoutput_PB3D() result(ierr)
            character(*), parameter :: rout_name = 'adapt_inoutput_PB3D'
            
            ! initialize ierr
            ierr = 0
            
            if (n_sol_requested.lt.1) then
                n_sol_requested = 1
                call writo('n_sol_requested has been increased to '&
                    &//trim(i2str(n_sol_requested)),warning=.true.)
            end if
            if (rich_restart_lvl.lt.1) then
                call writo('rich_restart_lvl was '//&
                    &trim(i2str(rich_restart_lvl))//&
                    &' but cannot be lower than 1, so it was reset to 1',&
                    &warning=.true.)
                rich_restart_lvl = 1
            end if
            if (rich_restart_lvl.gt.max_it_rich) then
                ierr = 1
                err_msg = 'rich_restart_lvl not within 1..max_it_rich = 1..'//&
                    &trim(i2str(max_it_rich))
                CHCKERR(err_msg)
            end if
            if (rich_restart_lvl.gt.1) then
                ierr = find_max_rich_lvl(PB3D_rich_lvl,PB3D_rich_lvl)
                CHCKERR('')
                if (rich_restart_lvl.gt.PB3D_rich_lvl+1) then
                    ierr = 1
                    if (PB3D_rich_lvl.gt.0) then
                        err_msg = 'The highest Richardson level found was '//&
                            &trim(i2str(PB3D_rich_lvl))
                    else
                        err_msg = 'No Richardson level found'
                    end if
                    CHCKERR(err_msg)
                end if
            end if
        end function adapt_inoutput_PB3D
        
        ! Checks whether the variables concerning plotting are chosen correctly:
        !   n_theta_plot and n_zeta_plot have to be positive,
        !   if n_theta_plot or  n_zeta_plot is equal to 1, the  minimum value is
        !   chosen.
        !   if  comparing different  toroidal positions,  only 2  values can  be
        !   used.
        !   ex_plot_style should be 1 (GNUPlot) or 2 (Bokeh / Mayavi)
        !   if toroidal positions are compared, theta_plot limits should contain
        !   0..2pi and max_r_plot should be 1.
        !> \private
        integer function adapt_plot() result(ierr)
            use num_vars, only: compare_tor_pos
            
            character(*), parameter :: rout_name = 'adapt_plot'
            
            ! local variables
            real(dp), parameter :: tol = 1.E-6_dp                               ! tolerance for checks
            character(len=max_str_ln) :: err_msg                                ! error message
            
            ! initialize ierr
            ierr = 0
            
            if (n_theta_plot.lt.1) then
                n_theta_plot = 1
                call writo('n_theta_plot has to be positive and is set to '//&
                    &trim(i2str(n_theta_plot)),warning=.true.)
            end if
            if (n_zeta_plot.lt.1) then
                n_zeta_plot = 1
                call writo('n_zeta_plot has to be positive and is set to '//&
                    &trim(i2str(n_zeta_plot)),warning=.true.)
            end if
            if (n_theta_plot.eq.1 .and. &
                &abs(max_theta_plot-min_theta_plot).gt.tol) then
                call writo('For n_theta_plot = 1, max_theta_plot is changed to &
                    &min_theta_plot = '//trim(r2strt(min_theta_plot)),&
                    &warning=.true.)
                max_theta_plot = min_theta_plot
            end if
            if (n_zeta_plot.eq.1 .and. &
                &abs(max_zeta_plot-min_zeta_plot).gt.tol) then
                call writo('For n_zeta_plot = 1, min_zeta_plot and &
                    &max_zeta_plot are changed to their average = '//&
                    &trim(r2strt((min_zeta_plot+max_zeta_plot)/2)),&
                    &warning=.true.)
                max_zeta_plot = 0.5*(min_zeta_plot+max_zeta_plot)
                min_zeta_plot = max_zeta_plot
            end if
            if (ex_plot_style.lt.1 .or. ex_plot_style.gt.2) then
                call writo('external plot style should be')
                call lvl_ud(1)
                call writo('1: GNUPlot')
                call writo('2: Bokeh (2-D) / Mayavi (3-D)')
                call lvl_ud(-1)
                call writo('Default (1: GNUPlot) chosen',warning=.true.)
                ex_plot_style = 1
            end if
            if (prog_style.eq.2) then                                           ! only for POST
                if (compare_tor_pos) then
                    if (n_zeta_plot.ne.3) then
                        ierr = 1
                        err_msg = 'For compare_tor_pos, you can only use &
                            &n_zeta = 3 (one point in the middle)'
                        CHCKERR(err_msg)
                    end if
                    if (POST_style.ne.1) then
                        ierr = 1
                        err_msg = 'For compare_tor_pos, you can only use &
                            &POST_style = 1 (extended grid)'
                        CHCKERR(err_msg)
                    end if
                    if (min_r_plot.lt.tol_zero) then
                        ierr = 1
                        err_msg = 'min_r_plot has to be greater than 0.'
                        CHCKERR(err_msg)
                    end if
                    if (max_r_plot.lt.1._dp) then
                        ierr = 1
                        err_msg = 'max_r_plot has to be 1.'
                        CHCKERR(err_msg)
                    end if
                    if (abs(min_theta_plot-0).gt.tol_zero .or. &
                        &abs(max_theta_plot-2).gt.tol_zero) then
                        ierr = 1
                        err_msg = 'theta limits of plot need to be 0..2pi'
                        CHCKERR(err_msg)
                    end if
                end if
                if (min_Rvac_plot.lt.0) then
                    call writo('The minimum value for Rvac_plot has to be &
                        &greater than zero',warning=.true.)
                    call lvl_ud(1)
                        min_Rvac_plot = 0.1*R_0
                        call writo('It has been set to '//&
                            &trim(r2str(min_Rvac_plot)))
                    call lvl_ud(-1)
                end if
                if (max_Rvac_plot.lt.min_Rvac_plot) then
                    ierr = 1
                    err_msg = 'The maximum value for Rvac_plot has to be &
                        &greater than the minimum value'
                    CHCKERR(err_msg)
                end if
                if (eq_style.eq.1 .and. plot_vac_pot) then
                    plot_vac_pot = .false.
                    call writo('Not plotting vacuum for VMEC yet',&
                        &warning=.true.)
                end if
            end if
        end function adapt_plot
        
        ! Checks whether  variables concerning  perturbation modes  are correct:
        !   the absolute value of prim_X has to be at least min_nm_X,
        !   deduces the X  style,
        !   min_sec_X  and  max_sec_X, pertaining  to  X  style 1  and  n_sec_X,
        !   pertaining to X style 2 cannot be both given,
        !   if  X  style  1  (prescribed),  max_sec_X has  to  be  greater  than
        !   min_sec_X in absolute value,
        !   if  X style 2 (fast), n_mod_X has to be a number greater than 0,
        !   for every X style, some checks are made and some variables set.
        !> \private
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
                    if (max_sec_X.ge.huge(1)) then
                        ierr = 1
                        err_msg = 'max_sec_X not set'
                        CHCKERR(err_msg)
                    end if
                    if (min_sec_X.ge.huge(1)) then
                        ierr = 1
                        err_msg = 'min_sec_X not set'
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
            end select
        end function adapt_X_modes
        
        ! Checks whether  min_n_par_X has  been chosen  correctly:
        !   It should be larger than  req_min_n_par_X (arbitrary) if an integral
        !   is to be calculated. However,  there are additional requirements for
        !   the different magnetic integral styles: for style 1 min_par_X has to
        !   be of the form 2+k and for style 2 of the form 4+3k.
        !> \private
        subroutine adapt_min_n_par_X
            ! local variables
            integer :: fund_n_par                                               ! fundamental interval width
            
            ! check req_min_n_par_X
            if (min_n_par_X.lt.req_min_n_par_X)  then
                min_n_par_X = req_min_n_par_X
                call writo('min_n_par_X has been increased to '//&
                    &trim(i2str(min_n_par_X))//' to have enough parallel &
                    &resolution',warning=.true.)
            end if
            
            ! check for individual magnetic integral style
            select case (magn_int_style)
                case (1)                                                        ! Trapezoidal rule
                    fund_n_par = 1
                case (2)                                                        ! Simpson's 3/8 rule
                    fund_n_par = 3
            end select
            if (mod(min_n_par_X-1,fund_n_par).ne.0) then                        ! there is a remainder
                min_n_par_X = 1 + fund_n_par * &
                    &((min_n_par_X-1)/fund_n_par + 1)
                call writo('min_n_par_X has been increased to '//&
                    &trim(i2str(min_n_par_X))//' to be a of the form '//&
                    &trim(i2str(fund_n_par+1))//' + '//&
                    &trim(i2str(fund_n_par))//&
                    &'k for magnetic integral style '//&
                    &trim(i2str(magn_int_style)),warning=.true.)
            end if
        end subroutine adapt_min_n_par_X
        
        ! Checks whether variables related to alpha are chosen correctly :
        !   alpha_style has to be either 0 or 1
        !   for 1, n_alpha = 1 and min/max_alpha = alpha
        !   for 2, max_par_X - min_par_X = 2
        !> \private
        integer function adapt_alpha() result(ierr)
            character(*), parameter :: rout_name = 'adapt_alpha'
            
            ! local variables
            character(len=max_str_ln) :: err_msg                                ! error message
            
            ! initialize ierr
            ierr = 0
            
            ! check alpha style
            if (alpha_style.lt.1 .or. alpha_style.gt.2) then
                ierr = 1
                err_msg = 'Alpha style needs to be 1 or 2 [def]'
                CHCKERR(err_msg)
            end if
            
            ! adapt variables pertaining to particular alpha style
            select case (alpha_style)
                case (1)                                                        ! single field line, multiple turns
                    if (n_alpha.ne.1) then
                        call writo('For alpha style 1, n_alpha can only be 1',&
                            &warning=.true.)
                        n_alpha = 1
                    end if
                    max_alpha = alpha
                case (2)                                                        ! multiple field lines, single turns
                    if (n_alpha.lt.1) then
                        ierr = 1
                        err_msg = 'For alpha style 2, n_alpha needs to be &
                            &positive'
                        CHCKERR(err_msg)
                    end if
                    if (abs(max_par_X - min_par_X - 2) .gt. tol_zero) then
                        min_par_X = 0._dp
                        max_par_X = 2._dp
                        call writo('min/max_par_X has been set to ['//&
                            &trim(r2strt(min_par_X))//':'//&
                            &trim(r2strt(max_par_X))//']',warning=.true.)
                    end if
            end select
            
            ! for HELENA, warn about axisymmetry
            if (eq_style.eq.2) then
                if (abs(max_par_X - min_par_X - 2) .gt. tol_zero) then
                    call writo('HELENA equilibria are axisymmetric, so only &
                        &max - min_par_X = 2 needed',warning=.true.)
                end if
                if (n_alpha .gt. 1) then
                    call writo('HELENA equilibria are axisymmetric, so only &
                        &n_alpha = 1 needed',warning=.true.)
                end if
            end if
        end function adapt_alpha
        
        ! Checks  whether variables  concerning the  solution grid  are correct:
        !   min_r should be greater than 0
        !   max_r should be lesser than 1
        !   min_r should be lesser than max_r
        ! The number of solution points should be large enough for the numerical
        ! discretization scheme, and there should be enough points per process.
        !> \private
        integer function adapt_sol_grid(min_r,max_r,var_name) result(ierr)
            use num_vars, only: n_procs
            
            character(*), parameter :: rout_name = 'adapt_sol_grid'
            
            ! input / output
            real(dp), intent(inout) :: min_r, max_r                             ! min and max of norma range
            character(len=*), intent(in) :: var_name                            ! name of variable
            
            ! local variables
            character(len=max_str_ln) :: err_msg                                ! error message
            
            ! initialize ierr
            ierr = 0
            
            ! check min_r
            if (min_r.lt.0._dp) then
                min_r = 0._dp
                call writo('min_r_'//trim(var_name)//' has been increased to '&
                    &//trim(r2strt(min_r)),warning=.true.)
            end if
            
            ! check if max_r is not greater than 1
            if (max_r.gt.1._dp) then
                max_r = 1._dp
                call writo('max_r_'//trim(var_name)//' has been decreased to '&
                    &//trim(r2strt(max_r)),warning=.true.)
            end if
            
            ! check if min_r < max_r
            if (min_r.ge.max_r) then
                ierr = 1
                err_msg = 'min_r_'//trim(var_name)//' should be lower than &
                    & max_r_'//trim(var_name)
                CHCKERR(err_msg)
            end if
            
            ! check n_r_sol
            if (n_r_sol.lt.6*norm_disc_prec_sol+2) then
                n_r_sol = 6*norm_disc_prec_sol+2
                call writo('n_r_sol has been increased to '//&
                    &trim(i2str(n_r_sol))//' to have enough points for the &
                    &normal discretization precision',warning=.true.)
            end if
            if (n_r_sol.lt.n_procs*(2*norm_disc_prec_sol+1)) then
                n_r_sol = n_procs*(2*norm_disc_prec_sol+1)
                call writo('n_r_sol has been increased to '//&
                    &trim(i2str(n_r_sol))//' to have enough points per &
                    &process',warning=.true.)
            end if
            
            ! warning if too close to axis for VMEC
            if (eq_style.eq.1) then
                if (min_r.le.tol_zero) &
                    &call writo('For VMEC, you should not go all the way to &
                    &the magnetic axis.',warning=.true.)
            end if
        end function adapt_sol_grid
        
        ! Checks  whether the variables concerning  Richardson extrapolation are
        ! correct:
        !   max_it_rich has to be at least 1.
        !> \private
        subroutine adapt_rich
            if (max_it_rich.lt.1) then
                max_it_rich = 1
                call writo('max_it_rich has been increased to 1',warning=.true.)
            end if
        end subroutine adapt_rich
        
        ! Checks whether the variables concerning finding zeros are correct:
        !   max_it_zero has to be at least 2,
        !   relax_fac_HH has to be larger than 0.
        !> \private
        subroutine adapt_zero
            if (max_it_zero.lt.1) then
                max_it_zero = 2
                call writo('max_it_zero has been increased to 2',warning=.true.)
            end if
            if (relax_fac_HH.lt.0) then
                relax_fac_HH = def_relax_fac_HH
                call writo('reset relax_fac_HH to '//&
                    &trim(r2strt(def_relax_fac_HH))&
                    &//' as it should be larger than 0',warning=.true.)
            end if
        end subroutine adapt_zero
        
        ! Checks whether tolerances are correct:
        !   tol_norm needs to be 0..1,
        !   tol_zero needs to be min_tol..max_tol,
        !   tol_rich needs to be min_tol..max_tol,
        !   tol_SLEPC needs to be min_tol..max_tol.
        ! Also sets local tol_SLEPC.
        !> \private
        subroutine adapt_tols
            ! local variables
            integer :: id                                                       ! counter
            integer :: max_id                                                   ! maximum of user-set id
            real(dp), parameter :: tol_SLEPC_def = 1.E-5_dp                     ! default tol_SLEPC
            
            ! check tol_norm
            if (tol_norm.lt.0._dp) then
                call writo('tol_norm has been increased to 0',warning=.true.)
                tol_norm = 0._dp
            end if
            if (tol_norm.gt.1._dp) then
                call writo('tol_norm has been decreased to 1',warning=.true.)
                tol_norm = 1._dp
            end if
            
            ! check tol_rich
            if (tol_rich.lt.min_tol) then
                call writo('tol_rich has been increased to '//&
                    &trim(r2str(min_tol)),warning=.true.)
                tol_rich = min_tol
            end if
            if (tol_rich.gt.max_tol) then
                call writo('tol_rich has been decreased to '//&
                    &trim(r2str(max_tol)),warning=.true.)
                tol_rich = max_tol
            end if
            
            ! check tol_zero
            if (tol_zero.lt.min_tol) then
                call writo('tol_zero has been increased to '//&
                    &trim(r2str(min_tol)),warning=.true.)
                tol_zero = min_tol
            end if
            if (tol_zero.gt.max_tol) then
                call writo('tol_zero has been decreased to '//&
                    &trim(r2str(max_tol)),warning=.true.)
                tol_zero = max_tol
            end if
            
            ! check and setup tol_SLEPC
            allocate(tol_SLEPC_loc(max_it_rich))
            ! get maximum index set by user and correct
            max_id = 0
            do id = 1,max_size_tol_SLEPC
                if (tol_SLEPC(id).lt.huge(1._dp)) then                          ! was overwritten by user
                    max_id = id
                    if (tol_SLEPC(id).lt.min_tol) then
                        call writo('tol_SLEPC('//trim(i2str(id))//&
                            &') has been increased to '//trim(r2str(min_tol)),&
                            &warning=.true.)
                        tol_SLEPC(id) = min_tol
                    end if
                    if (tol_SLEPC(id).gt.max_tol) then
                        call writo('tol_SLEPC('//trim(i2str(id))//&
                            &') has been decreased to '//trim(r2str(max_tol)),&
                            &warning=.true.)
                        tol_SLEPC(id) = max_tol
                    end if
                end if
            end do
            ! set local variable
            if (max_id.gt.0) then                                               ! was overwritten by user
                if (max_id.gt.max_it_rich) then                                 ! too many values given
                    call writo(trim(i2str(max_id-max_it_rich))//&
                        &' last values discarded for tol_SLEPC',warning=.true.)
                    tol_SLEPC_loc = tol_SLEPC(1:max_it_rich)
                else if (max_id.eq.max_it_rich) then                            ! exactly the right amount of values given
                    tol_SLEPC_loc = tol_SLEPC(1:max_it_rich)
                else                                                            ! not enough values given
                    call writo('tol_SLEPC('//trim(i2str(max_id+1))//':'//&
                        &trim(i2str(max_it_rich))//') has been set to '//&
                        &trim(r2str(tol_SLEPC_def))//' [default]',&
                        &warning=.true.)
                    tol_SLEPC_loc(1:max_id) = tol_SLEPC(1:max_id)
                    tol_SLEPC_loc(max_id+1:max_it_rich) = tol_SLEPC_def
                end if
            else                                                                ! was not overwritten by user
                tol_SLEPC_loc = tol_SLEPC_def                                   ! all tolerances equal to default
            end if
        end subroutine adapt_tols
        
        ! Checks whether normalization variables are chosen correctly:
        !   rho_0 has to be positive.
        !> \private
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
    
    !> Reads the equilibrium input file if no Richardson restart.
    !!
    !! These variables will be passed on through the HDF5 data system.
    !!
    !! \see reconstruct_pb3d_in().
    !!
    !! \return ierr
    integer function read_input_eq() result(ierr)
        use num_vars, only: eq_style, use_pol_flux_E, eq_i
        use VMEC_ops, only: read_VMEC
        use HELENA_ops, only: read_HEL
        use grid_vars, only: n_r_in
        
        character(*), parameter :: rout_name = 'read_input_eq'
        
        
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
        end select
        
        ! close equilibrium file
        close(eq_i)
    end function read_input_eq
    
    !> Print input quantities to an output file.
    !!
    !! Variables printed:
    !!  -  \c misc_in:  \c  prog_version,  \c eq_style,  \c  rho_style, \c  R_0,
    !!  \c  pres_0,   \c  B_0,  \c  psi_0   \c  rho_0,  \c  T_0,   \c  vac_perm,
    !!  \c   use_pol_flux_F,  \c   use_pol_flux_E,   \c  use_normalization,   \c
    !!  norm_disc_prec_eq, \c norm_disc_style_sol, \c n_r_in, \c n_r_eq, \c
    !!  n_r_sol, \c debug_version
    !!  -  \c misc_in_V: \c is_asym_V,  \c is_freeb_V, \c mnmax_V,  \c mpol_V, \c
    !!  ntor_V, \c gam_V
    !!  - \c flux_t_V
    !!  - \c flux_p_V
    !!  - \c pres_V
    !!  - \c rot_t_V
    !!  - \c q_saf_V
    !!  - \c mn_V
    !!  - \c R_V_c
    !!  - \c R_V_s
    !!  - \c Z_V_c
    !!  - \c Z_V_s
    !!  - \c L_V_c
    !!  - \c L_V_s
    !!  - \c B_V_sub: \c B_V_sub_c, \c B_V_sub_s
    !!  - \c B_V: \c B_V_c, \c B_V_s
    !!  - \c J_V_sup_int
    !!  - \c jac_V: \c jac_V_c, \c jac_V_s
    !!  - \c misc_in_H: \c ias, \c nchi
    !!  - \c RZ_H: \c R_H, \c Z_H
    !!  - \c chi_H
    !!  - \c q_saf_H
    !!  - \c rot_t_H
    !!  - \c RBphi_H
    !!  - \c pres_H
    !!  - \c flux_p_H
    !!  - \c flux_t_H
    !!  -  \c  misc_X:  \c  prim_X,  \c n_mod_X,  \c  min_sec_X,  \c  max_sec_X,
    !!  \c  norm_disc_prec_X,   \c  norm_style,  \c  U_style,   \c  X_style,  \c
    !!  matrix_SLEPC_style, \c K_style, \c alpha_style
    !!  - \c  misc_sol: \c  min_r_sol, \c  max_r_sol, \c  norm_disc_prec_sol, \c
    !!  BC_style, \c EV_BC, \c EV_BC
    !!
    !! \return ierr
    integer function print_output_in(data_name) result(ierr)
        use num_vars, only: eq_style, rho_style, prog_version, use_pol_flux_E, &
            &use_pol_flux_F, use_normalization, norm_disc_prec_eq, PB3D_name, &
            &norm_disc_prec_X, norm_style, U_style, X_style, alpha_style, &
            &matrix_SLEPC_style, BC_style, EV_style, norm_disc_prec_sol, &
            &EV_BC, magn_int_style, K_style, debug_version, plot_VMEC_modes, &
            &norm_disc_style_sol, X_grid_style, V_interp_style
        use eq_vars, only: R_0, pres_0, B_0, psi_0, rho_0, T_0, vac_perm, &
            &max_flux_E, max_flux_F
        use grid_vars, onLy: n_r_in, n_r_eq, n_r_sol, n_alpha, min_alpha, &
            &max_alpha
        use grid_ops, only: calc_norm_range
        use X_vars, only: min_r_sol, max_r_sol, min_sec_X, max_sec_X, prim_X, &
            &n_mod_X
        use HDF5_ops, only: print_HDF5_arrs
        use HDF5_vars, only: dealloc_var_1D, var_1D_type, &
            &max_dim_var_1D
        use HELENA_vars, only: chi_H, flux_p_H, flux_t_H, R_H, Z_H, nchi, ias, &
            &q_saf_H, rot_t_H, pres_H, RBphi_H
        use VMEC_vars, only: is_freeb_V, mnmax_V, mpol_V, ntor_V, is_asym_V, &
            &gam_V, R_V_c, R_V_s, Z_V_c, Z_V_s, L_V_c, L_V_s, jac_V_c, &
            &jac_V_s, mnmax_V, mn_V, rot_t_V, q_saf_V, pres_V, flux_t_V, &
            &flux_p_V, nfp_V
        use HELENA_vars, only: h_H_11, h_H_12, h_H_33
#if ldebug
        use VMEC_vars, only: B_V_sub_c, B_V_sub_s, B_V_c, B_V_s, J_V_sup_int
#endif
        
        character(*), parameter :: rout_name = 'print_output_in'
        
        ! input / output
        character(len=*), intent(in) :: data_name                               !< name under which to store
        
        ! local variables
        type(var_1D_type), allocatable, target :: in_1D(:)                      ! 1D equivalent of input variables
        type(var_1D_type), pointer :: in_1D_loc => null()                       ! local element in in_1D
        integer :: id                                                           ! counter
        integer :: in_limits(2)                                                 ! min. and max. index of input variable grid of this process
        integer :: loc_size                                                     ! local size
        
        ! initialize ierr
        ierr = 0
        
        ! user output
        call writo('Write input variables to output file')
        call lvl_ud(1)
        
        ! plot modes if requested and VMEC
        if (plot_VMEC_modes) then
            call writo('Plotting decay of VMEC modes')
            call print_ex_2D('R_V_c','R_V_c',&
                &log10(maxval(abs(R_V_c(:,:,0)),2)),draw=.false.)
            call draw_ex(['R_V_c'],'R_V_c',1,1,.false.,ex_plot_style=1)
            call print_ex_2D('R_V_s','R_V_s',&
                &log10(maxval(abs(R_V_s(:,:,0)),2)),draw=.false.)
            call draw_ex(['R_V_s'],'R_V_s',1,1,.false.,ex_plot_style=1)
            call print_ex_2D('Z_V_c','Z_V_c',&
                &log10(maxval(abs(Z_V_c(:,:,0)),2)),draw=.false.)
            call draw_ex(['Z_V_c'],'Z_V_c',1,1,.false.,ex_plot_style=1)
            call print_ex_2D('Z_V_s','Z_V_s',&
                &log10(maxval(abs(Z_V_s(:,:,0)),2)),draw=.false.)
            call draw_ex(['Z_V_s'],'Z_V_s',1,1,.false.,ex_plot_style=1)
            call print_ex_2D('L_V_c','L_V_c',&
                &log10(maxval(abs(L_V_c(:,:,0)),2)),draw=.false.)
            call draw_ex(['L_V_c'],'L_V_c',1,1,.false.,ex_plot_style=1)
            call print_ex_2D('L_V_s','L_V_s',&
                &log10(maxval(abs(L_V_s(:,:,0)),2)),draw=.false.)
            call draw_ex(['L_V_s'],'L_V_s',1,1,.false.,ex_plot_style=1)
            call print_ex_2D('jac_V_c','jac_V_c',&
                &log10(maxval(abs(jac_V_c(:,:,0)),2)),draw=.false.)
            call draw_ex(['jac_V_c'],'jac_V_c',1,1,.false.,ex_plot_style=1)
            call print_ex_2D('jac_V_s','jac_V_s',&
                &log10(maxval(abs(jac_V_s(:,:,0)),2)),draw=.false.)
            call draw_ex(['jac_V_s'],'jac_V_s',1,1,.false.,ex_plot_style=1)
#if ldebug
            call print_ex_2D('B_V_c','B_V_c',&
                &log10(maxval(abs(B_V_c),2)),draw=.false.)
            call draw_ex(['B_V_c'],'B_V_c',1,1,.false.,ex_plot_style=1)
            call print_ex_2D('B_V_s','B_V_s',&
                &log10(maxval(abs(B_V_s),2)),draw=.false.)
            call draw_ex(['B_V_s'],'B_V_s',1,1,.false.,ex_plot_style=1)
#endif
        end if
        
        ! calculate limits of input range
        ierr = calc_norm_range(in_limits=in_limits)
        CHCKERR('')
        n_r_eq = in_limits(2)-in_limits(1)+1
        if (in_limits(1).gt.1 .or. in_limits(2).lt.n_r_in) then
            call writo('Only the normal range '//trim(i2str(in_limits(1)))//&
                &'..'//trim(i2str(in_limits(2)))//' is written')
        else
            call writo('The entire normal range '//trim(i2str(in_limits(1)))//&
                &'..'//trim(i2str(in_limits(2)))//' is written')
        end if
        call writo('(relevant to solution range '//&
            &trim(r2strt(min_r_sol))//'..'//trim(r2strt(max_r_sol))//')')
        
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
        in_1D_loc%loc_i_max = [20]
        in_1D_loc%tot_i_min = in_1D_loc%loc_i_min
        in_1D_loc%tot_i_max = in_1D_loc%loc_i_max
        allocate(in_1D_loc%p(20))
        in_1D_loc%p = [prog_version,eq_style*1._dp,rho_style*1._dp,&
            &R_0,pres_0,B_0,psi_0,rho_0,T_0,vac_perm,-1._dp,-1._dp,&
            &-1._dp,norm_disc_prec_eq*1._dp,n_r_in*1._dp,n_r_eq*1._dp,&
            &n_r_sol*1._dp,max_flux_E,max_flux_F,-1._dp]
        if (use_pol_flux_E) in_1D_loc%p(11) = 1._dp
        if (use_pol_flux_F) in_1D_loc%p(12) = 1._dp
        if (use_normalization) in_1D_loc%p(13) = 1._dp
        if (debug_version) in_1D_loc%p(20) = 1._dp
        
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
                
                loc_size = n_r_eq*size(flux_t_V,2)
                
                ! flux_t_V
                in_1D_loc => in_1D(id); id = id+1
                in_1D_loc%var_name = 'flux_t_V'
                allocate(in_1D_loc%tot_i_min(2),in_1D_loc%tot_i_max(2))
                allocate(in_1D_loc%loc_i_min(2),in_1D_loc%loc_i_max(2))
                in_1D_loc%loc_i_min = [1,0]
                in_1D_loc%loc_i_max = [n_r_eq,size(flux_t_V,2)-1]
                in_1D_loc%tot_i_min = in_1D_loc%loc_i_min
                in_1D_loc%tot_i_max = in_1D_loc%loc_i_max
                allocate(in_1D_loc%p(loc_size))
                in_1D_loc%p = reshape(flux_t_V(in_limits(1):in_limits(2),:),&
                    &[loc_size])
                
                ! flux_p_V
                in_1D_loc => in_1D(id); id = id+1
                in_1D_loc%var_name = 'flux_p_V'
                allocate(in_1D_loc%tot_i_min(2),in_1D_loc%tot_i_max(2))
                allocate(in_1D_loc%loc_i_min(2),in_1D_loc%loc_i_max(2))
                in_1D_loc%loc_i_min = [1,0]
                in_1D_loc%loc_i_max = [n_r_eq,size(flux_p_V,2)-1]
                in_1D_loc%tot_i_min = in_1D_loc%loc_i_min
                in_1D_loc%tot_i_max = in_1D_loc%loc_i_max
                allocate(in_1D_loc%p(loc_size))
                in_1D_loc%p = reshape(flux_p_V(in_limits(1):in_limits(2),:),&
                    &[loc_size])
                
                loc_size = n_r_eq*size(pres_V,2)
                
                ! pres_V
                in_1D_loc => in_1D(id); id = id+1
                in_1D_loc%var_name = 'pres_V'
                allocate(in_1D_loc%tot_i_min(2),in_1D_loc%tot_i_max(2))
                allocate(in_1D_loc%loc_i_min(2),in_1D_loc%loc_i_max(2))
                in_1D_loc%loc_i_min = [1,0]
                in_1D_loc%loc_i_max = [n_r_eq,size(pres_V,2)-1]
                in_1D_loc%tot_i_min = in_1D_loc%loc_i_min
                in_1D_loc%tot_i_max = in_1D_loc%loc_i_max
                allocate(in_1D_loc%p(loc_size))
                in_1D_loc%p = reshape(pres_V(in_limits(1):in_limits(2),:),&
                    &[loc_size])
                
                ! rot_t_V
                in_1D_loc => in_1D(id); id = id+1
                in_1D_loc%var_name = 'rot_t_V'
                allocate(in_1D_loc%tot_i_min(2),in_1D_loc%tot_i_max(2))
                allocate(in_1D_loc%loc_i_min(2),in_1D_loc%loc_i_max(2))
                in_1D_loc%loc_i_min = [1,0]
                in_1D_loc%loc_i_max = [n_r_eq,size(rot_t_V,2)-1]
                in_1D_loc%tot_i_min = in_1D_loc%loc_i_min
                in_1D_loc%tot_i_max = in_1D_loc%loc_i_max
                allocate(in_1D_loc%p(loc_size))
                in_1D_loc%p = reshape(rot_t_V(in_limits(1):in_limits(2),:),&
                    &[loc_size])
                
                ! q_saf_V
                in_1D_loc => in_1D(id); id = id+1
                in_1D_loc%var_name = 'q_saf_V'
                allocate(in_1D_loc%tot_i_min(2),in_1D_loc%tot_i_max(2))
                allocate(in_1D_loc%loc_i_min(2),in_1D_loc%loc_i_max(2))
                in_1D_loc%loc_i_min = [1,0]
                in_1D_loc%loc_i_max = [n_r_eq,size(q_saf_V,2)-1]
                in_1D_loc%tot_i_min = in_1D_loc%loc_i_min
                in_1D_loc%tot_i_max = in_1D_loc%loc_i_max
                allocate(in_1D_loc%p(loc_size))
                in_1D_loc%p = reshape(q_saf_V(in_limits(1):in_limits(2),:),&
                    &[loc_size])
                
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
                allocate(in_1D_loc%tot_i_min(4),in_1D_loc%tot_i_max(4))
                allocate(in_1D_loc%loc_i_min(4),in_1D_loc%loc_i_max(4))
                in_1D_loc%loc_i_min = [1,1,0,1]
                in_1D_loc%loc_i_max = [mnmax_V,n_r_eq,size(R_V_c,3)-1,6]
                in_1D_loc%tot_i_min = in_1D_loc%loc_i_min
                in_1D_loc%tot_i_max = in_1D_loc%loc_i_max
                allocate(in_1D_loc%p(6*mnmax_V*n_r_eq*size(R_V_c,3)))
                in_1D_loc%p = reshape([R_V_c(:,in_limits(1):in_limits(2),:),&
                    &R_V_s(:,in_limits(1):in_limits(2),:),&
                    &Z_V_c(:,in_limits(1):in_limits(2),:),&
                    &Z_V_s(:,in_limits(1):in_limits(2),:),&
                    &L_V_c(:,in_limits(1):in_limits(2),:),&
                    &L_V_s(:,in_limits(1):in_limits(2),:)],&
                    &[6*mnmax_V*n_r_eq*size(R_V_c,3)]) 
                
                ! jac_V
                in_1D_loc => in_1D(id); id = id+1
                in_1D_loc%var_name = 'jac_V'
                allocate(in_1D_loc%tot_i_min(4),in_1D_loc%tot_i_max(4))
                allocate(in_1D_loc%loc_i_min(4),in_1D_loc%loc_i_max(4))
                in_1D_loc%loc_i_min = [1,1,0,1]
                in_1D_loc%loc_i_max = [mnmax_V,n_r_eq,size(jac_V_c,3)-1,2]
                in_1D_loc%tot_i_min = in_1D_loc%loc_i_min
                in_1D_loc%tot_i_max = in_1D_loc%loc_i_max
                allocate(in_1D_loc%p(2*mnmax_V*n_r_eq*size(jac_V_c,3)))
                in_1D_loc%p = reshape([jac_V_c(:,in_limits(1):in_limits(2),:),&
                    &jac_V_s(:,in_limits(1):in_limits(2),:)],&
                    &[2*mnmax_V*n_r_eq*size(jac_V_c,3)])
                
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
                
                ! J_V_sup_int
                in_1D_loc => in_1D(id); id = id+1
                in_1D_loc%var_name = 'J_V_sup_int'
                allocate(in_1D_loc%tot_i_min(2),in_1D_loc%tot_i_max(2))
                allocate(in_1D_loc%loc_i_min(2),in_1D_loc%loc_i_max(2))
                in_1D_loc%loc_i_min = [1,1]
                in_1D_loc%loc_i_max = [n_r_eq,2]
                in_1D_loc%tot_i_min = in_1D_loc%loc_i_min
                in_1D_loc%tot_i_max = in_1D_loc%loc_i_max
                allocate(in_1D_loc%p(2*n_r_eq))
                in_1D_loc%p = reshape(J_V_sup_int(in_limits(1):in_limits(2),:),&
                    &[2*n_r_eq])
                
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
                
                loc_size = n_r_eq*size(flux_p_H,2)
                
                ! flux_p_H
                in_1D_loc => in_1D(id); id = id+1
                in_1D_loc%var_name = 'flux_p_H'
                allocate(in_1D_loc%tot_i_min(2),in_1D_loc%tot_i_max(2))
                allocate(in_1D_loc%loc_i_min(2),in_1D_loc%loc_i_max(2))
                in_1D_loc%loc_i_min = [1,0]
                in_1D_loc%loc_i_max = [n_r_eq,size(flux_p_H,2)-1]
                in_1D_loc%tot_i_min = in_1D_loc%loc_i_min
                in_1D_loc%tot_i_max = in_1D_loc%loc_i_max
                allocate(in_1D_loc%p(loc_size))
                in_1D_loc%p = reshape(flux_p_H(in_limits(1):in_limits(2),:),&
                    &[loc_size])
                
                ! flux_t_H
                in_1D_loc => in_1D(id); id = id+1
                in_1D_loc%var_name = 'flux_t_H'
                allocate(in_1D_loc%tot_i_min(2),in_1D_loc%tot_i_max(2))
                allocate(in_1D_loc%loc_i_min(2),in_1D_loc%loc_i_max(2))
                in_1D_loc%loc_i_min = [1,0]
                in_1D_loc%loc_i_max = [n_r_eq,size(flux_t_H,2)-1]
                in_1D_loc%tot_i_min = in_1D_loc%loc_i_min
                in_1D_loc%tot_i_max = in_1D_loc%loc_i_max
                allocate(in_1D_loc%p(loc_size))
                in_1D_loc%p = reshape(flux_t_H(in_limits(1):in_limits(2),:),&
                    &[loc_size])
                
                ! q_saf_H
                in_1D_loc => in_1D(id); id = id+1
                in_1D_loc%var_name = 'q_saf_H'
                allocate(in_1D_loc%tot_i_min(2),in_1D_loc%tot_i_max(2))
                allocate(in_1D_loc%loc_i_min(2),in_1D_loc%loc_i_max(2))
                in_1D_loc%loc_i_min = [1,0]
                in_1D_loc%loc_i_max = [n_r_eq,size(q_saf_H,2)-1]
                in_1D_loc%tot_i_min = in_1D_loc%loc_i_min
                in_1D_loc%tot_i_max = in_1D_loc%loc_i_max
                allocate(in_1D_loc%p(loc_size))
                in_1D_loc%p = reshape(q_saf_H(in_limits(1):in_limits(2),:),&
                    &[loc_size])
                
                ! rot_t_H
                in_1D_loc => in_1D(id); id = id+1
                in_1D_loc%var_name = 'rot_t_H'
                allocate(in_1D_loc%tot_i_min(2),in_1D_loc%tot_i_max(2))
                allocate(in_1D_loc%loc_i_min(2),in_1D_loc%loc_i_max(2))
                in_1D_loc%loc_i_min = [1,0]
                in_1D_loc%loc_i_max = [n_r_eq,size(rot_t_H,2)-1]
                in_1D_loc%tot_i_min = in_1D_loc%loc_i_min
                in_1D_loc%tot_i_max = in_1D_loc%loc_i_max
                allocate(in_1D_loc%p(loc_size))
                in_1D_loc%p = reshape(rot_t_H(in_limits(1):in_limits(2),:),&
                    &[loc_size])
                
                ! pres_H
                in_1D_loc => in_1D(id); id = id+1
                in_1D_loc%var_name = 'pres_H'
                allocate(in_1D_loc%tot_i_min(2),in_1D_loc%tot_i_max(2))
                allocate(in_1D_loc%loc_i_min(2),in_1D_loc%loc_i_max(2))
                in_1D_loc%loc_i_min = [1,0]
                in_1D_loc%loc_i_max = [n_r_eq,size(pres_H,2)-1]
                in_1D_loc%tot_i_min = in_1D_loc%loc_i_min
                in_1D_loc%tot_i_max = in_1D_loc%loc_i_max
                allocate(in_1D_loc%p(loc_size))
                in_1D_loc%p = reshape(pres_H(in_limits(1):in_limits(2),:),&
                    &[loc_size])
                
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
                
                ! h_H
                in_1D_loc => in_1D(id); id = id+1
                in_1D_loc%var_name = 'h_H'
                allocate(in_1D_loc%tot_i_min(3),in_1D_loc%tot_i_max(3))
                allocate(in_1D_loc%loc_i_min(3),in_1D_loc%loc_i_max(3))
                in_1D_loc%loc_i_min = [1,1,1]
                in_1D_loc%loc_i_max = [nchi,n_r_eq,4]
                in_1D_loc%tot_i_min = in_1D_loc%loc_i_min
                in_1D_loc%tot_i_max = in_1D_loc%loc_i_max
                allocate(in_1D_loc%p(3*nchi*n_r_eq))
                in_1D_loc%p = reshape([h_H_11(:,in_limits(1):in_limits(2)),&
                    &h_H_12(:,in_limits(1):in_limits(2)),&
                    &h_H_33(:,in_limits(1):in_limits(2))],&
                    &[3*nchi*n_r_eq])
        end select
        
        ! misc_X
        in_1D_loc => in_1D(id); id = id+1
        in_1D_loc%var_name = 'misc_X'
        allocate(in_1D_loc%tot_i_min(1),in_1D_loc%tot_i_max(1))
        allocate(in_1D_loc%loc_i_min(1),in_1D_loc%loc_i_max(1))
        in_1D_loc%loc_i_min = [1]
        in_1D_loc%loc_i_max = [17]
        in_1D_loc%tot_i_min = in_1D_loc%loc_i_min
        in_1D_loc%tot_i_max = in_1D_loc%loc_i_max
        allocate(in_1D_loc%p(17))
        in_1D_loc%p = [prim_X*1._dp,n_mod_X*1._dp,min_sec_X*1._dp,&
            &max_sec_X*1._dp,norm_disc_prec_X*1._dp,norm_style*1._dp,&
            &U_style*1._dp,X_style*1._dp,matrix_SLEPC_style*1._dp,&
            &magn_int_style*1._dp,K_style*1._dp,alpha_style*1._dp,&
            &X_grid_style*1._dp,V_interp_style*1._dp,min_alpha*1._dp,&
            &max_alpha*1._dp,n_alpha*1._dp]
        
        ! misc_sol
        in_1D_loc => in_1D(id); id = id+1
        in_1D_loc%var_name = 'misc_sol'
        allocate(in_1D_loc%tot_i_min(1),in_1D_loc%tot_i_max(1))
        allocate(in_1D_loc%loc_i_min(1),in_1D_loc%loc_i_max(1))
        in_1D_loc%loc_i_min = [1]
        in_1D_loc%loc_i_max = [8]
        in_1D_loc%tot_i_min = in_1D_loc%loc_i_min
        in_1D_loc%tot_i_max = in_1D_loc%loc_i_max
        allocate(in_1D_loc%p(8))
        in_1D_loc%p = [min_r_sol,max_r_sol,norm_disc_prec_sol*1._dp,&
            &norm_disc_style_sol*1._dp,BC_style(1)*1._dp,BC_style(2)*1._dp,&
            &EV_style*1._dp,EV_BC]
        
        ! write
        ierr = print_HDF5_arrs(in_1D(1:id-1),PB3D_name,trim(data_name))
        CHCKERR('')
        
        ! clean up
        call dealloc_var_1D(in_1D)
        nullify(in_1D_loc)
        
        ! user output
        call lvl_ud(-1)
    end function print_output_in
end module input_ops
