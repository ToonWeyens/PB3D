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
    public read_input, pause_prog, get_real, get_int, get_log, test_max_memory
    
contains
    ! Queries for a logical value yes or no, where the default answer is also to
    ! be provided.
    ! [MPI] All ranks, but only global rank can give input
    function get_log(yes,ind) result(val)
        use messages, only: start_time, stop_time
        use num_vars, only: glb_rank
        use MPI_utilities, only: wait_MPI, broadcast_var
        
        ! input / output
        logical :: val                                                          ! output value
        logical :: yes
        logical, intent(in), optional :: ind                                    ! individual pause or not
        
        ! local variables
        character(len=11) :: empty_str = ''                                     ! empty string
        character(len=max_str_ln) :: answer_str                                 ! string with answer
        integer :: istat                                                        ! status
        logical :: ind_loc                                                      ! local version of ind
        
        ! set local ind
        ind_loc = .false.
        if (present(ind)) ind_loc = ind
        
        ! only global master can receive input
        if (glb_rank.eq.0) then
            write(*,'(A)',advance='no') empty_str                               ! first print empty string so that output is visible
            if (yes) then
                write(*,'(A)',advance='no') 'y(es)/n(o) [yes]: '
                val = .true.
            else
                write(*,'(A)',advance='no') 'y(es)/n(o) [no]: '
                val = .false.
            end if
            call stop_time
            read (*, '(A)') answer_str
            call start_time
            
            select case (strh2l(trim(answer_str)))
                !case ('y','yes')
                case ('y','Y','yes','Yes','YEs','YES','yEs','yES','yeS','YeS')
                    val = .true.
                case ('n','N','no','No','NO','nO') 
                    val = .false.
                case default 
            end select
        end if
        
        ! if not individual, broadcast result
        if (.not.ind_loc) then
            istat = broadcast_var(val)
            if (istat.ne.0) call writo('WARNING: In get_log, something went &
                &wrong. Default used.')
        end if
    end function get_log
    
    ! Queries  for user  input for a  real value, where  allowable range  can be
    ! provided as well.
    ! [MPI] All ranks, but only global rank can give input
    function get_real(lim_lo,lim_hi,ind) result(val)
        use messages, only: start_time, stop_time
        use num_vars, only: glb_rank
        use MPI_utilities, only: wait_MPI, broadcast_var
        
        ! input / output
        real(dp) :: val                                                         ! output value
        real(dp), intent(in), optional :: lim_lo                                ! upper and lower limit
        real(dp), intent(in), optional :: lim_hi                                ! upper and lower limit
        logical, intent(in), optional :: ind                                    ! individual pause or not
        
        ! local variables
        character(len=11) :: empty_str = ''                                     ! empty string
        integer :: istat                                                        ! status
        logical :: ind_loc                                                      ! local version of ind
        real(dp) :: lims_loc(2)                                                 ! local version of limits
        
        ! set local ind
        ind_loc = .false.
        if (present(ind)) ind_loc = ind
        
        ! initialize val
        val = 0._dp
        
        ! only global master can receive input
        if (glb_rank.eq.0) then
            ! set up local limits
            lims_loc = [-huge(1._dp),huge(1._dp)]
            if (present(lim_lo)) lims_loc(1) = lim_lo
            if (present(lim_hi)) lims_loc(2) = lim_hi
            
            ! get input
            do
                write(*,'(A)',advance='no') empty_str                               ! first print empty string so that output is visible
                write(*,'(A)',advance='no') 'Input a value'
                if (present(lim_lo).or.present(lim_hi)) then
                    write(*,'(A)',advance='no') ' ['
                    if (present(lim_lo)) write(*,'(A)',advance='no') &
                        &trim(r2strt(lim_lo))
                    write(*,'(A)',advance='no') '..'
                    if (present(lim_hi)) write(*,'(A)',advance='no') &
                        &trim(r2strt(lim_hi))
                    write(*,'(A)',advance='no') ']'
                end if
                write(*,'(A)',advance='no') ': '
                call stop_time
                read (*,*,iostat=istat) val
                call start_time
                
                if (istat.ne.0 .or. val.lt.lims_loc(1) .or. val.gt. &
                    &lims_loc(2)) then
                    write(*,'(A)',advance='no') empty_str                       ! first print empty string so that output is visible
                    write(*,'(A)',advance='no') 'Choose a value between '//&
                        &trim(r2strt(lims_loc(1)))//' and '//&
                        &trim(r2strt(lims_loc(2)))
                    write(*,*) ''
                    cycle
                else
                    exit
                end if
            end do
        end if
        
        ! if not individual, broadcast result
        if (.not.ind_loc) then
            istat = broadcast_var(val)
            if (istat.ne.0) call writo('WARNING: In get_real, something went &
                &wrong. Default of zero used.')
        end if
    end function get_real
    
    ! Queries for user input for an  integer value, where allowable range can be
    ! provided as well.
    ! [MPI] All ranks, but only global rank can give input
    function get_int(lim_lo,lim_hi,ind) result(val)
        use messages, only: start_time, stop_time
        use num_vars, only: glb_rank
        use MPI_utilities, only: wait_MPI, broadcast_var
        
        ! input / output
        integer :: val                                                          ! output value
        integer, intent(in), optional :: lim_lo                                 ! upper and lower limit
        integer, intent(in), optional :: lim_hi                                 ! upper and lower limit
        logical, intent(in), optional :: ind                                    ! individual pause or not
        
        ! local variables
        character(len=11) :: empty_str = ''                                     ! empty string
        integer :: istat                                                        ! status
        logical :: ind_loc                                                      ! local version of ind
        integer :: lims_loc(2)                                                ! local version of limits
        
        ! set local ind
        ind_loc = .false.
        if (present(ind)) ind_loc = ind
        
        ! initialize val
        val = 0
        
        ! only global master can receive input
        if (glb_rank.eq.0) then
            ! set up local limits
            lims_loc = [-huge(1),huge(1)]
            if (present(lim_lo)) lims_loc(1) = lim_lo
            if (present(lim_hi)) lims_loc(2) = lim_hi
            
            ! get input
            do
                write(*,'(A)',advance='no') empty_str                               ! first print empty string so that output is visible
                write(*,'(A)',advance='no') 'Input a value'
                if (present(lim_lo).or.present(lim_hi)) then
                    write(*,'(A)',advance='no') ' ['
                    if (present(lim_lo)) write(*,'(A)',advance='no') &
                        &trim(i2str(lim_lo))
                    write(*,'(A)',advance='no') '..'
                    if (present(lim_hi)) write(*,'(A)',advance='no') &
                        &trim(i2str(lim_hi))
                    write(*,'(A)',advance='no') ']'
                end if
                write(*,'(A)',advance='no') ': '
                call stop_time
                read (*,*,iostat=istat) val
                call start_time
                
                if (istat.ne.0 .or. val.lt.lims_loc(1) .or. val.gt. &
                    &lims_loc(2)) then
                    write(*,'(A)',advance='no') empty_str                       ! first print empty string so that output is visible
                    write(*,'(A)',advance='no') 'Choose a value between '//&
                        &trim(i2str(lims_loc(1)))//' and '//&
                        &trim(i2str(lims_loc(2)))
                    write(*,*) ''
                    cycle
                else
                    exit
                end if
            end do
        end if
        
        ! if not individual, broadcast result
        if (.not.ind_loc) then
            istat = broadcast_var(val)
            if (istat.ne.0) call writo('WARNING: In get_int, something went &
                &wrong. Default of zero used.')
        end if
    end function get_int
    
    ! pauses the running of the program
    ! [MPI] All ranks or, optinally, only current rank
    subroutine pause_prog(ind)
        use messages, only: start_time, stop_time
        use MPI_utilities, only: wait_MPI
        use num_vars, only: glb_rank, grp_rank
        
        ! input / output
        logical, intent(in), optional :: ind                                    ! individual pause or not
        
        ! local variables
        character(len=11) :: empty_str = ''                                     ! empty string
        integer :: istat                                                        ! status
        logical :: ind_loc                                                      ! local version of ind
        
        ! output message
        if (grp_rank.eq.0) then
            write(*,'(A)',advance='no') empty_str                               ! first print empty string so that output is visible
            write(*,'(A)',advance='no') 'Paused. Press enter...'
        end if
        
        ! set local ind
        ind_loc = .false.
        if (present(ind)) ind_loc = ind
        
        ! only global master can receive input
        if (glb_rank.eq.0) then
            call stop_time
            read (*, *)
            call start_time
        end if
        
        ! wait for MPI
        if (.not.ind_loc) then
            istat = wait_MPI()
            if (istat.ne.0) call writo('WARNING: In pause_prog, something went &
                &wrong. Continuing.')
        end if
    end subroutine pause_prog
    
    ! reads input from user-provided input file
    ! [MPI] only global master
    integer function read_input() result(ierr)
        use num_vars, only: &
            &minim_style, min_alpha, max_alpha, n_alpha, max_it_NR, tol_NR, &
            &max_it_r, input_i, use_pol_flux_F, EV_style, max_mem_per_proc, &
            &n_procs_per_alpha, plot_resonance, tol_r, n_sol_requested, &
            &nyq_fac, glb_rank, nyq_fac, plot_grid, plot_flux_q, &
            &use_normalization, n_sol_plotted, n_theta_plot, n_zeta_plot, &
            &EV_BC, rho_style, retain_all_sol, prog_style, norm_disc_ord, &
            &BC_style, max_it_inv, tol_norm_r, tol_slepc, max_n_it_slepc, &
            &glb_n_procs
        use eq_vars, only: rho_0
        use messages, only: writo, lvl_ud
        use files_ops, only: input_name
        use X_vars, only: min_n_X, max_n_X, min_m_X, max_m_X, min_n_r_X, &
            &min_n_X, min_r_X, max_r_X
        use grid_vars, only: n_par_X, min_par_X, max_par_X
        
        character(*), parameter :: rout_name = 'read_input'
        
        ! local variables
        integer :: istat                                                        ! error
        integer :: prim_X, min_sec_X, max_sec_X                                 ! n_X and m_X (pol. flux) or m_X and n_X (tor. flux)
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! input options
        namelist /inputdata_PB3D/ max_mem_per_proc, minim_style, min_par_X, &
            &min_n_r_X, max_par_X, min_alpha, max_alpha, n_par_X, n_alpha, &
            &max_it_NR, tol_NR, max_it_r, tol_r, prim_X, min_sec_X, max_sec_X, &
            &min_r_X, max_r_X, EV_style, n_procs_per_alpha, plot_resonance, &
            &n_sol_requested, EV_BC, tol_slepc, rho_style, nyq_fac, rho_0, &
            &use_pol_flux_F, plot_grid, plot_flux_q, use_normalization, &
            &n_theta_plot, n_zeta_plot, retain_all_sol, norm_disc_ord, &
            &BC_style, max_it_inv, tol_norm_r, max_n_it_slepc
        namelist /inputdata_PB3D_POST/ max_mem_per_proc, n_sol_plotted, &
            &n_theta_plot, n_zeta_plot, plot_resonance, plot_flux_q, plot_grid
        
        ! initialize ierr
        ierr = 0
        
        if (glb_rank.eq.0) then                                                 ! only global master
            if (input_name.ne.'') then                                          ! if open_input opened a file
                call writo('Setting up user-provided input "' &       
                    &// trim(input_name) // '"')
            else 
                call writo('Setting default values for input')
            end if
            call lvl_ud(1)
            
            ! initialize input variables (optionally overwritten by user later)
            if (input_name.ne.'') call writo('Initialize all the inputs with &
                &default values')
            
            ! common variables for all program styles
            max_mem_per_proc = 6000_dp/glb_n_procs                              ! count with 6GB
            
            ! select depending on program style
            select case (prog_style)
                case(1)                                                         ! PB3D
                    call default_input_PB3D
                case(2)                                                         ! PB3D_POST
                    call default_input_PB3D_POST
                case default
                    err_msg = 'No program style associated with '//&
                        &trim(i2str(prog_style))
                    ierr = 1
                    CHCKERR(err_msg)
            end select
            
            ! read user input
            if (input_name.ne.'') then                                          ! otherwise, defaults are loaded
                ! select depending on program style
                select case (prog_style)
                    case(1)                                                     ! PB3D
                        read(input_i,nml=inputdata_PB3D,iostat=istat)           ! read input data
                        
                        ! check input if successful read
                        if (istat.eq.0) then                                    ! input file succesfully read
                            call writo('Overwriting with user-provided file "'&
                                &//trim(input_name) // '"')
                            
                            call writo('Checking user-provided file')
                            
                            call lvl_ud(1)
                            
                            ! adapt run-time variables if needed
                            call adapt_run
                            
                            ! adapt output variables if needed
                            call adapt_output
                            
                            ! adapt plotting variables if needed
                            call adapt_plot
                            
                            ! adapt alpha variables if needed
                            call adapt_alpha
                            
                            ! adapt n_par_X if needed
                            call adapt_n_par_X
                            
                            ! adapt min_r_X and max_r_X if needed
                            ierr = adapt_X()
                            CHCKERR('')
                            
                            ! adapt Richardson variables if needed
                            call adapt_r
                            
                            ! adapt variables for inverse if needed
                            call adapt_inv
                            
                            ! adapt Newton-Rhapson variables if needed
                            call adapt_NR
                            
                            ! adapt m variables if needed
                            ierr = adapt_m()
                            CHCKERR('')
                            
                            ! adapt normalization variables if needed
                            ierr = adapt_normalization()
                            CHCKERR('')
                            
                            call lvl_ud(-1)
                        else                                                    ! cannot read input data
                            call writo('WARNING: Cannot open user-provided &
                                &file "'//trim(input_name)//'". Using defaults')
                        end if
                        
                        ! set up min_n_X, max_n_X, min_m_X, max_m_X
                        if (use_pol_flux_F) then
                            min_n_X = prim_X
                            max_n_X = prim_X
                            min_m_X = min_sec_X
                            max_m_X = max_sec_X
                        else
                            min_m_X = prim_X
                            max_m_X = prim_X
                            min_n_X = min_sec_X
                            max_n_X = max_sec_X
                        end if
                    case(2)                                                     ! PB3D_POST
                        read(input_i,nml=inputdata_PB3D_POST,iostat=istat)      ! read input data
                        
                        ! check input if successful read
                        if (istat.eq.0) then                                    ! input file succesfully read
                            call writo('Overwriting with user-provided file "'&
                                &//trim(input_name) // '"')
                        else                                                    ! cannot read input data
                            call writo('WARNING: Cannot open user-provided &
                                &file "'//trim(input_name)//'". Using defaults')
                        end if
                    case default
                        err_msg = 'No program style associated with '//&
                            &trim(i2str(prog_style))
                        ierr = 1
                        CHCKERR(err_msg)
                end select
                
                ! user output
                call writo('Close input file')
                close(input_i)
            end if
            
            call lvl_ud(-1)
            call writo('Input values set')
        end if
    contains
        subroutine default_input_PB3D
            use num_vars, only: eq_style, use_pol_flux_E
            
            ! concerning Newton-Rhapson
            max_it_NR = 500                                                     ! maximum 500 Newton-Rhapson iterations
            tol_NR = 1.0E-10_dp                                                 ! wanted relative error in Newton-Rhapson iteration
            
            ! concerning Richardson extrapolation
            max_it_r = 1                                                        ! by default no Richardson extrapolation
            tol_r = 1E-5                                                        ! wanted relative error in Richardson extrapolation
            
            ! concerning calculating the inverse
            max_it_inv = 1                                                      ! by default no iteration to calculate inverse
            
            ! runtime variables
            minim_style = 1                                                     ! Richardson Extrapolation with normal discretization
            n_procs_per_alpha = 1                                               ! 1 processor per field line
            max_n_it_slepc = 1000                                               ! max. nr. of iterations for SLEPC
            plot_resonance = .false.                                            ! do not plot the q-profile with nq-m = 0
            plot_grid = .false.                                                 ! do not plot the grid
            plot_flux_q = .false.                                               ! do not plot the flux quantities
            use_normalization = .true.                                          ! use normalization for the variables
            EV_BC = 1._dp                                                       ! use 1 as artificial EV for the Boundary Conditions
            tol_slepc = 1.E-8_dp                                                ! tolerance of 1E-8
            rho_style = 1                                                       ! constant pressure profile, equal to rho_0
            norm_disc_ord = 1                                                   ! order 1 normal discretization
            BC_style = [1,2]                                                    ! left BC zeroed and right BC through minimization of energy
            
            ! variables concerning output
            n_sol_requested = 3                                                 ! request solutions with 3 highes EV
            retain_all_sol = .false.                                            ! don't retain faulty ones
            ! default   values   of  n_theta_plot  and   n_zeta_plot  depend  on
            ! equilibrium style being used:
            !   1:  VMEC
            !   2:  HELENA
            select case(eq_style)
                case (1)                                                        ! VMEC
                    n_theta_plot = 201                                          ! nr. poloidal points in plot
                    n_zeta_plot = 101                                           ! nr. toroidal points in plot
                case (2)                                                        ! HELENA
                    n_theta_plot = 501                                          ! nr. poloidal points in plot
                    n_zeta_plot = 1                                             ! nr. toroidal points in plot
            end select
            
            ! variables concerning poloidal mode numbers m
            min_par_X = -4.0_dp                                                 ! minimum parallel angle [pi]
            max_par_X = 4.0_dp                                                  ! maximum parallel angle [pi]
            nyq_fac = 10                                                        ! need at least 10 points per period for perturbation quantitites
            prim_X = 20                                                         ! main mode number of perturbation
            min_sec_X = prim_X                                                  ! min. of. secondary mode number of perturbation
            max_sec_X = prim_X                                                  ! max. of. secondary mode number of perturbation
            n_par_X = 20                                                        ! number of parallel grid points in pert. grid
            use_pol_flux_F = use_pol_flux_E                                     ! use same normal flux coordinate as the equilibrium
            
            ! variables concerning alpha
            min_alpha = 0.0_dp                                                  ! minimum field line label [pi]
            max_alpha = 2.0_dp                                                  ! maximum field line label [pi]
            n_alpha = 10                                                        ! number of different field lines
            tol_norm_r = 0.05                                                   ! tolerance for normal range
            
            ! variables concerning perturbation
            min_r_X = 0.1_dp                                                    ! minimum radius
            max_r_X = 1.0_dp                                                    ! maximum radius
            EV_style = 1                                                        ! slepc solver for EV problem
            min_n_r_X = 20                                                      ! at least 20 points in perturbation grid
            
            ! variables concerning normalization
            rho_0 = 10E-6_dp                                                    ! for fusion, particle density of around 1E21, mp around 1E-27
        end subroutine default_input_PB3D
        
        subroutine default_input_PB3D_POST
            use num_vars, only: eq_style
            
            ! runtime variables
            plot_resonance = .false.                                            ! do not plot the q-profile with nq-m = 0
            plot_flux_q = .false.                                               ! do not plot the flux quantities
            plot_grid = .false.                                                 ! do not plot the grid
            
            ! variables concerning output
            n_sol_plotted = n_sol_requested                                     ! plot all solutions
            ! default   values   of  n_theta_plot  and   n_zeta_plot  depend  on
            ! equilibrium style being used:
            !   1:  VMEC
            !   2:  HELENA
            select case(eq_style)
                case (1)                                                        ! VMEC
                    n_theta_plot = 201                                          ! nr. poloidal points in plot
                    n_zeta_plot = 101                                           ! nr. toroidal points in plot
                case (2)                                                        ! HELENA
                    n_theta_plot = 501                                          ! nr. poloidal points in plot
                    n_zeta_plot = 1                                             ! nr. toroidal points in plot
            end select
        end subroutine default_input_PB3D_POST
        
        ! checks whether the variables concerning run-time are chosen correctly.
        ! n_procs_per_alpha has to be at least 1
        ! rho_style has to be 1 (constant rho = rho_0)
        subroutine adapt_run
            if (n_procs_per_alpha.lt.1) then
                n_procs_per_alpha = 1
                call writo('WARNING: n_procs_per_alpha has been increased to '&
                    &//trim(i2str(n_procs_per_alpha)))
            end if
            if (rho_style.ne.1) then
                rho_style = 1
                call writo('WARNING: rho_style set to default (1: constant)')
            end if
        end subroutine adapt_run
        
        ! checks whether the variables concerning output are chosen correctly.
        ! n_sol_requested has to be at least one
        subroutine adapt_output
            if (n_sol_requested.lt.1) then
                n_sol_requested = 1
                call writo('WARNING: n_sol_requested has been increased to '&
                    &//trim(i2str(n_sol_requested)))
            end if
        end subroutine adapt_output
        
        ! checks whether the variables concerning plotting are chosen correctly.
        ! n_theta and n_zeta_plot have to be positive
        subroutine adapt_plot
            if (n_theta_plot.lt.1) then
                n_theta_plot = 1
                call writo('WARNING: n_theta_plot cannot be negative and is &
                    &set to '//trim(i2str(n_theta_plot)))
            end if
            if (n_zeta_plot.lt.1) then
                n_zeta_plot = 1
                call writo('WARNING: n_zeta_plot cannot be negative and is &
                    &set to '//trim(i2str(n_zeta_plot)))
            end if
        end subroutine adapt_plot
        
        ! checks  whether n_par_X  is  chosen  high enough  so  aliasing can  be
        ! avoided.  aliasing occurs  when there  are  not enough  points on  the
        ! parallel grid so the fast-moving functions e^(i(k-m)) V don't give the
        ! wrong integrals in the perturbation part
        subroutine adapt_n_par_X
            if (n_par_X.lt.nyq_fac*max(max_sec_X-min_sec_X,1)*&
                &(max_par_X-min_par_X)/2) then
                n_par_X = int(nyq_fac*max(max_sec_X-min_sec_X,1)*&
                    &(max_par_X-min_par_X)/2)
                if (mod(n_par_X,2).eq.0) n_par_X = n_par_X+1                    ! odd numbers are usually better
                call writo('WARNING: To avoid aliasing of the perturbation &
                    &integrals, n_par_X is increased to '//trim(i2str(n_par_X)))
            end if
        end subroutine adapt_n_par_X
        
        ! checks whether variables concerning  perturbation are correct. min_r_X
        ! should not be  too close to zero because  the equilibrium calculations
        ! yield an infinity at the magnetic  axis. max_r_X cannot be larger than
        ! 1.0 and has to  be larger than min_r_X
        ! the absolute  value of prim_X  has to be  at least 5  (arbitrary), but
        ! preferibly at least 10 (arbitrary)
        ! min_n_r_X has to be at  least 6*norm_disc_ord+2 (for at least two full
        ! points)
        integer function adapt_X() result(ierr)
            use grid_vars, only: n_r_eq
            
            character(*), parameter :: rout_name = 'adapt_X'
            
            ! local variables
            character(len=max_str_ln) :: err_msg                                ! error message
            real(dp) :: one = 1.00000001_dp                                     ! one plus a little margin
            integer :: abs_min_prim_X = 5                                       ! absolute minimum for the high-n theory
            integer :: rec_min_prim_X = 10                                      ! recomended minimum for the high-n theory
            
            ! initialize ierr
            ierr = 0
            
            ! check prim_X
            if (abs(prim_X).lt.abs_min_prim_X) then
                prim_X = sign(abs_min_prim_X,prim_X)
                call writo('WARNING: prim_X has been changed to '//&
                    &trim(i2str(abs_min_prim_X))//' but its absolute value &
                    &should be at least '//trim(i2str(rec_min_prim_X))//&
                    &' for correct results')
            else if (abs(prim_X).lt.rec_min_prim_X) then
                call writo('WARNING: The absolute value of prim_X should be at &
                    &least '//trim(i2str(rec_min_prim_X))//' for correct &
                    &results')
            end if
            
            ! check min_n_r_X
            if (min_n_r_X.lt.6*norm_disc_ord+2) then
                min_n_r_X = 6*norm_disc_ord+2
                call writo('WARNING: min_n_r_X has been increased to '//&
                    &trim(i2str(min_n_r_X)))
            end if
            
            ! check min_r_X
            if (min_r_X.lt.one/(n_r_eq-1)) then
                min_r_X = one/(n_r_eq-1)
                call writo('WARNING: min_r_X has been increased to '//&
                    &trim(r2strt(min_r_X)))
            end if
            
            ! check if max_r_X is not greater than 1
            if (max_r_X.gt.1.0) then
                max_r_X = 1.
                call writo('WARNING: max_r_X has been decreased to '//&
                    &trim(r2strt(max_r_X)))
            end if
            
            ! check  if min_r_X  < max_r_X  with at  least one  equilbrium point
            ! between them
            if (min_r_X+1./(n_r_eq-1).ge.max_r_X) then
                ierr = 1
                err_msg = 'max_r_X - min_r_X has to be at least '//&
                    &trim(r2strt(1._dp/(n_r_eq-1)))
                CHCKERR(err_msg)
            end if
        end function adapt_X
        
        ! checks  whether the variables concerning  Richardson extrapolation are
        ! correct. max_it_r has to be at least 1
        subroutine adapt_r
            if (max_it_r.lt.1) then
                max_it_r = 1
                call writo('WARNING: max_it_r has been increased to 1')
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
        
        ! checks whether the variables concerning the poloidal mode number m are
        ! correct. max_sec_X has to be greater than min_sec_X and nyq_fac has to
        ! be at least 1
        ! sets up min_n_X, max_n_X, min_m_X, max_m_X using prim_X and min_sec_X,
        ! max_sec_X
        integer function adapt_m() result (ierr)
            character(*), parameter :: rout_name = 'adapt_m'
            
            ! local variables
            character(len=max_str_ln) :: err_msg                                ! error message
            
            ! initialize ierr
            ierr = 0
            
            ! check min_sec_X, max_sec_X
            if (max_sec_X.lt.min_sec_X) then
                ierr = 1
                err_msg = 'max_sec_X has to be larger or equal to min_sec_X'
                CHCKERR(err_msg)
            end if
            
            ! adapt nyq_fac
            if (nyq_fac.lt.1) then
                call writo('WARNING: nyq_fac has been increased to 1')
                nyq_fac = 1
            end if
        end function adapt_m
        
        ! checks  whether  n_alpha  is  chosen high  enough  and whether  normal
        ! tolerance is within 0 and 1
        subroutine adapt_alpha
            if (n_alpha.lt.1) then
                call writo('WARNING: n_alpha has been increased to 1')
                n_alpha = 1
            end if
            if (tol_norm_r.lt.0) then
                call writo('WARNING: tol_norm_r has been increased to 0')
                tol_norm_r = 0._dp
            end if
            if (tol_norm_r.gt.1) then
                call writo('WARNING: tol_norm_r has been decreased to 1')
                tol_norm_r = 1._dp
            end if
        end subroutine adapt_alpha
        
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
    end function read_input

    ! test whether maximum memory feasible
    integer function test_max_memory() result(ierr)
        use num_vars, only: max_mem_per_proc
        
        character(*), parameter :: rout_name = 'test_max_memory'
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer :: n_max                                                        ! maximum size of array
        real(dp), allocatable :: max_mem_arr(:,:)                               ! array with maximum size
        
        ! initialize ierr
        ierr = 0
        
        call writo('Testing whether maximum memory per process of '//&
            &trim(r2strt(max_mem_per_proc/1000))//'GB is possible')
        
        call lvl_ud(1)
        
        ! (lazy) allocation
        n_max = ceiling(sqrt(max_mem_per_proc/(sizeof(1._dp)*1.E-6)))
        call writo('Allocating doubles array of size '//trim(i2str(n_max))&
            &//'x'//trim(i2str(n_max)))
        allocate(max_mem_arr(n_max,n_max),STAT=ierr)
        err_msg = 'cannot allocate this much memory. Try setting &
            &"max_mem_per_proc" lower'
        CHCKERR(err_msg)
        
        ! explicitely set elements
        max_mem_arr = 0._dp                                                     ! this can fail while lazy allocation does not
        
        deallocate(max_mem_arr)
        
        call lvl_ud(-1)
        call writo('Maximum memory allocatable')
        
    end function test_max_memory
end module input_ops
