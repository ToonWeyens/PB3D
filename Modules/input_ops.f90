!------------------------------------------------------------------------------!
!   This module contains operations concerning giving input                    !
!------------------------------------------------------------------------------!
module input_ops
#include <PB3D_macros.h>
    use str_ops, only: strh2l, i2str, r2strt
    use num_vars, only: dp, max_str_ln
    
    implicit none
    private
    public yes_no, read_input

contains
    ! queries for yes or no, depending on the flag yes:
    !   yes = .true.: yes is default answer
    !   yes = .false.: no is default answer
    logical function yes_no(yes)
        use message_ops, only: start_time, stop_time, &
            &lvl_sep, lvl
        
        ! input / output
        logical :: yes
        
        ! local variables
        character(len=max_str_ln) :: answer_str
        integer :: id 

        do id = 1,lvl                                                           ! to get the response on the right collumn
            write(*,'(A)',advance='no') lvl_sep
        end do
        if (yes) then
            write(*,'(A)',advance='no') 'y(es)/n(o) [yes]: '
            yes_no = .true.
        else
            write(*,'(A)',advance='no') 'y(es)/n(o) [no]: '
            yes_no = .false.
        end if
        call stop_time
        read (*, '(A)') answer_str
        call start_time

        select case (strh2l(trim(answer_str)))
            !case ('y','Y','yes','Yes','YEs','YES','yEs','yES','yeS','YeS')
            case ('y','yes')
                yes_no = .true.
            case ('n','N','no','No','NO','nO') 
                yes_no = .false.
            case default 
        end select
    end function yes_no
    
    ! reads input from user-provided input file
    ! [MPI] only global master
    integer function read_input() result(ierr)
        use num_vars, only: &
            &minim_style, min_alpha, max_alpha, n_alpha, max_it_NR, tol_NR, &
            &max_it_r, input_i, n_seq_0, min_n_r_X, use_pol_flux, &
            &calc_mesh_style, EV_style, n_procs_per_alpha, plot_jq, tol_r, &
            &n_sol_requested, min_r_X, max_r_X, nyq_fac, max_n_plots, &
            &glb_rank, nyq_fac, plot_grid, output_style
        use eq_vars, only: &
            &min_par, max_par, n_par, rho_0
        use message_ops, only: writo, lvl_ud
        use file_ops, only: input_name
        use X_vars, only: min_n_X, max_n_X, min_m_X, max_m_X
        
        character(*), parameter :: rout_name = 'read_input'
        
        ! local variables
        integer :: istat                                                        ! error
        integer :: prim_X, min_sec_X, max_sec_X                                 ! n_X and m_X (pol. flux) or m_X and n_X (tor. flux)
        
        ! input options
        namelist /inputdata/ minim_style, min_par, min_n_r_X, &
            &max_par, min_alpha, max_alpha, n_par, n_alpha, max_it_NR, tol_NR, &
            &max_it_r, tol_r, prim_X, min_sec_X, max_sec_X, min_r_X, &
            &max_r_X, EV_style, n_procs_per_alpha, plot_jq, n_sol_requested, &
            &nyq_fac, rho_0, max_n_plots, use_pol_flux, plot_grid, output_style
        
        ! initialize ierr
        ierr = 0
        
        if (glb_rank.eq.0) then                                                 ! only global master
            if (input_i.ge.0) then                                              ! if open_input opened a file
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
            
            ! initialize non-input file variables
            calc_mesh_style = 0                                                 ! for debugging, not for users
            
            ! read user input
            if (input_i.ge.n_seq_0) then                                        ! otherwise, defaults are loaded
                read (input_i, nml=inputdata, iostat=istat)                     ! read input data
                
                ! check input if successful read
                if (istat.eq.0) then                                            ! input file succesfully read
                    call writo('Overwriting with user-provided file "' // &
                        &trim(input_name) // '"')
                    
                    call writo('Checking user-provided file')
                    
                    call lvl_ud(1)
                    
                    ! adapt run-time variables if needed
                    call adapt_run
                    
                    ! adapt alpha variables if needed
                    call adapt_n_alpha
                    
                    ! adapt n_par if needed
                    call adapt_n_par
                    
                    ! adapt min_r_X and max_r_X if needed
                    ierr = adapt_X()
                    CHCKERR('')
                    
                    ! adapt Richardson variables if needed
                    call adapt_r
                    
                    ! adapt Newton-Rhapson variables if needed
                    call adapt_NR
                    
                    ! adapt m variables if needed
                    ierr = adapt_m()
                    CHCKERR('')
                    
                    ! adapt normalization variables if needed
                    ierr = adapt_normalization()
                    CHCKERR('')
                    
                    call lvl_ud(-1)
                else                                                            ! cannot read input data
                    call writo('Cannot open user-provided file "' // &
                        &trim(input_name) // '". Using defaults')
                end if
            end if
            
            call lvl_ud(-1)
            call writo('Input values set')
        end if
    contains
        subroutine default_input
            use eq_vars, only: eq_use_pol_flux
            
            ! concerning Newton-Rhapson
            max_it_NR = 50                                                      ! maximum 50 Newton-Rhapson iterations
            tol_NR = 1.0E-10_dp                                                 ! wanted relative error in Newton-Rhapson iteration
            ! concerning Richardson extrapolation
            max_it_r = 8                                                        ! maximum 5 levels of Richardson extrapolation
            tol_r = 1E-5                                                        ! wanted relative error in Richardson extrapolation
            ! runtime variables
            minim_style = 1                                                     ! Richardson Extrapolation with normal discretization
            n_procs_per_alpha = 1                                               ! 1 processor per field line
            plot_jq = .false.                                                   ! do not plot the q-profile with nq-m = 0
            plot_grid = .false.                                                 ! do not plot the grid
            n_sol_requested = 3                                                 ! request solutions with 3 highes EV
            max_n_plots = 4                                                     ! maximum nr. of modes for which to plot output in plot_X_vec
            output_style = 1                                                    ! GNUPlot output
            ! variables concerning poloidal mode numbers m
            min_par = -4.0_dp                                                   ! minimum parallel angle [pi]
            max_par = 4.0_dp                                                    ! maximum parallel angle [pi]
            nyq_fac = 5                                                         ! need at least 5 points per period for perturbation quantitites
            prim_X = 20                                                         ! main mode number of perturbation
            min_sec_X = prim_X                                                  ! min. of. secondary mode number of perturbation
            max_sec_X = prim_X                                                  ! max. of. secondary mode number of perturbation
            n_par = 20                                                          ! number of parallel grid points
            use_pol_flux = eq_use_pol_flux                                      ! use same normal flux coordinate as the equilibrium
            ! variables concerning alpha
            min_alpha = 0.0_dp                                                  ! minimum field line label [pi]
            max_alpha = 2.0_dp                                                  ! maximum field line label [pi]
            n_alpha = 10                                                        ! number of different field lines
            ! variables concerning perturbation
            min_r_X = 0.1_dp                                                    ! minimum radius
            max_r_X = 1.0_dp                                                    ! maximum radius
            EV_style = 1                                                        ! slepc solver for EV problem
            min_n_r_X = 10                                                      ! at least 10 points in perturbation grid
            ! variables concerning normalization
            rho_0 = 10E-6_dp                                                    ! for fusion, particle density of around 1E21, mp around 1E-27
        end subroutine
        
        ! checks whether the variables concerning run-time are chosen correctly.
        ! n_procs_per_alpha and n_sol_requested have to be at least 1
        ! output_style has to be 1 (GNUPlot) or 2 (HDF5) (see output_ops)
        subroutine adapt_run
            if (n_procs_per_alpha.lt.1) then
                n_procs_per_alpha = 1
                call writo('WARNING: n_procs_per_alpha has been increased to '&
                    &//trim(i2str(n_procs_per_alpha)))
            end if
            if (n_sol_requested.lt.1) then
                n_sol_requested = 1
                call writo('WARNING: n_sol_requested has been increased to '&
                    &//trim(i2str(n_sol_requested)))
            end if
            if (max_n_plots.lt.0) then
                max_n_plots = 0
                call writo('WARNING: max_n_plots cannot be negative and is &
                    &set to 0')
            end if
            if (output_style.lt.1 .or. output_style.gt.2) then
                output_style = 1
                call writo('WARNING: output_style set to default (1: GNUPlot)')
            end if
        end subroutine adapt_run
        
        ! checks whether n_par is chosen high enough so aliasing can be avoided.
        ! aliasing occurs when there are not  enough points on the parallel grid
        ! so  the  fast-moving  functions  e^(i(k-m)) V  don't  give  the  wrong
        ! integrals in the perturbation part
        subroutine adapt_n_par
            use num_vars, only: eq_style
            use HEL_vars, only: nchi, ias
            
            ! local variables
            integer :: n_par_old                                                ! backup of n_par
            real(dp) :: min_par_old, max_par_old                                ! backup of min_par_old, max_par_old
            real(dp) :: tol = 1.0E-8_dp                                         ! tolerance for rounding to integers
            
            if (eq_style.eq.2) then                                             ! HELENA input: n_par is dictated by nchi
                ! set back ups
                n_par_old = n_par
                min_par_old = min_par
                max_par_old = max_par
                ! round min_par and max_par to nearest integers
                min_par = floor(min_par)
                max_par = ceiling(max_par)
                if (abs(max_par-min_par).lt.1._dp) max_par = max_par + 1
                if (abs(min_par-min_par_old).gt.tol .or. &
                    &abs(max_par-max_par_old).gt.tol) &
                    &call writo('WARNING: For HELENA output, min_par and &
                    &max_par are rounded to the integer values '//&
                    &trim(i2str(floor(min_par)))//' and '//&
                    &trim(i2str(ceiling(max_par))))
                if (ias.eq.0) then                                              ! top-bottom symmetry in HELENA
                    n_par = nint(max_par-min_par) * (nchi-1)                    ! nchi - 1 per half circle
                else                                                            ! no top-bottom symmetry in HELENA
                    n_par = nint(max_par-min_par) * nchi / 2                    ! nchi / 2 per half circle
                end if
                if (n_par.ne.n_par_old) then
                    call writo('WARNING: Nr. of parallel points dictated by &
                        &HELENA: '//trim(i2str(n_par)))
                end if
            else                                                                ! not HELENA input: n_par can be chosen independently
                if (n_par.lt.nyq_fac*(max_sec_X-min_sec_X)*&
                    &(max_par-min_par)/2) then
                    n_par = int(nyq_fac*(max_sec_X-min_sec_X)*&
                        &(max_par-min_par)/2)
                    call writo('WARNING: To avoid aliasing of the perturbation &
                        &integrals, n_par is increased to '//trim(i2str(n_par)))
                end if
            end if
        end subroutine adapt_n_par
        
        ! checks whether variables concerning  perturbation are correct. min_r_X
        ! should not be  too close to zero because  the equilibrium calculations
        ! yield an infinity at the magnetic  axis. max_r_X cannot be larger than
        ! 1.0 and has to  be larger than min_r_X
        ! the absolute  value of prim_X  has to be  at least 5  (arbitrary), but
        ! preferibly at least 10 (arbitrary)
        ! min_n_r_X has  to be at  least 2 (for 2 grid points)
        integer function adapt_X() result(ierr)
            use eq_vars, only: n_r_eq
            
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
            if (min_n_r_X.lt.2) then
                min_n_r_X = 2
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
            
            ! set up min_n_X, max_n_X, min_m_X, max_m_X
            if (use_pol_flux) then
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
        end function adapt_m
        
        ! checks whether n_alpha is chosen high enough
        subroutine adapt_n_alpha
            if (n_alpha.lt.1) then
                call writo('WARNING: n_alpha has been increased to 1')
                n_alpha = 1
            end if
        end subroutine adapt_n_alpha
        
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
end module input_ops
