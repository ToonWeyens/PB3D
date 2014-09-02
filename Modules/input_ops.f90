!------------------------------------------------------------------------------!
!   This module contains operations concerning giving input                    !
!------------------------------------------------------------------------------!
module input_ops
#include <PB3D_macros.h>
    use str_ops, only: strh2l, i2str, r2strt
    use num_vars, only: &
        &dp, max_str_ln, style, min_alpha, max_alpha, n_alpha, max_it_NR, &
        &tol_NR, max_it_r, theta_var_along_B, input_i, n_seq_0, min_n_r_X, &
        &calc_mesh_style, EV_style, n_procs_per_alpha, plot_q, tol_r, &
        &n_sol_requested, min_r_X, max_r_X, reuse_r, nyq_fac, pi
    use eq_vars, only: &
        &min_par, max_par, n_par
    use output_ops, only: writo, lvl_ud, &
        &format_out
    use file_ops, only: input_name
    use X_vars, only: n_X, min_m_X, max_m_X
    implicit none
    private
    public yes_no, read_input

    ! input options
    namelist /inputdata/ format_out, style, min_par, min_n_r_X, &
        &max_par, min_alpha, max_alpha, n_par, n_alpha, max_it_NR, tol_NR, &
        &max_it_r, tol_r, theta_var_along_B, n_X, min_m_X, max_m_X, min_r_X, &
        &max_r_X, EV_style, n_procs_per_alpha, plot_q, n_sol_requested, &
        &reuse_r, nyq_fac

contains
    ! queries for yes or no, depending on the flag yes:
    !   yes = .true.: yes is default answer
    !   yes = .false.: no is default answer
    logical function yes_no(yes)
        use output_ops, only: start_time, stop_time, &
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
        use num_vars, only: glb_rank, nyq_fac
        
        character(*), parameter :: rout_name = 'read_input'
        
        ! local variables
        integer :: istat                                                        ! error
        
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
                if (istat.eq.0) then                                            ! input file succesfully read
                    call writo('Overwriting with user-provided file "' // &
                        &trim(input_name) // '"')
                    
                    call writo('Checking user-provided file')
                    
                    call lvl_ud(1)
                    
                    ! adapt run-time variables if needed
                    call adapt_run
                    
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
                    
                    ! adapt alpha variables if needed
                    ierr = adapt_alpha()
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
            use num_vars, only: pi
            
            ! concerning Newton-Rhapson
            max_it_NR = 50                                                      ! maximum 50 Newton-Rhapson iterations
            tol_NR = 1.0E-10_dp                                                 ! wanted relative error in Newton-Rhapson iteration
            ! concerning Richardson extrapolation
            max_it_r = 8                                                        ! maximum 5 levels of Richardson extrapolation
            tol_r = 1E-5                                                        ! wanted relative error in Richardson extrapolation
            ! runtime variables
            reuse_r = .true.                                                    ! reuse A and B from previous Richardson level
            format_out = 1                                                      ! NETCDF output
            style = 1                                                           ! Richardson Extrapolation with normal discretization
            n_procs_per_alpha = 1                                               ! 1 processor per field line
            plot_q = .false.                                                    ! do not plot the q-profile with nq-m = 0
            n_sol_requested = 3                                                 ! request solutions with 3 highes EV
            ! variables concerning poloidal mode numbers m
            min_m_X = 20                                                        ! lowest poloidal mode number m_X
            max_m_X = 22                                                        ! highest poloidal mode number m_X
            min_par = -4.0_dp*pi                                                ! minimum parallel angle
            max_par = 4.0_dp*pi                                                 ! maximum parallel angle
            nyq_fac = 5                                                         ! need at least 5 points per period for perturbation quantitites
            n_par = nyq_fac*(max_m_X-min_m_X)                                   ! number of parallel grid points
            ! variables concerning alpha
            min_alpha = 0.0_dp                                                  ! minimum field line label
            max_alpha = 2.0_dp*pi                                               ! maximum field line label
            n_alpha = 10                                                        ! number of different field lines
            theta_var_along_B = .true.                                          ! theta is used as the default parallel variable
            ! variables concerning perturbation
            n_X = 20                                                            ! toroidal mode number of perturbation
            min_r_X = 0.1_dp                                                    ! minimum radius
            max_r_X = 1.0_dp                                                    ! maximum radius
            EV_style = 1                                                        ! slepc solver for EV problem
            min_n_r_X = 10                                                      ! at least 10 points in perturbation grid
        end subroutine
        
        ! checks whether the variables concerning run-time are chosen correctly.
        ! n_procs_per_alpha and n_sol_requested have to be at least 1
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
        end subroutine adapt_run
        
        ! checks whether n_par is chosen high enough so aliasing can be avoided.
        ! aliasing occurs when there are not  enough points on the parallel grid
        ! so  the  fast-moving  functions  e^(i(k-m)) V  don't  give  the  wrong
        ! integrals in the perturbation part
        subroutine adapt_n_par
            if (n_par.lt.nyq_fac*(max_m_X-min_m_X)) then
                n_par = nyq_fac*(max_m_X-min_m_X)
                call writo('WARNING: To avoid aliasing of the perturbation &
                    &integrals, n_par is increased to '//trim(i2str(n_par)))
            end if
        end subroutine adapt_n_par
        
        ! checks whether variables concerning  perturbation are correct. min_r_X
        ! should not be  too close to zero because  the equilibrium calculations
        ! yield an infinity at the magnetic  axis. max_r_X cannot be larger than
        ! 1.0  and has  to be  larger than  min_r_X. n_X  has to  be at  least 5
        ! (arbitrary), but preferibly at least  10 (arbitrary). min_n_r_X has to
        ! be at least 2 (for 2 grid points)
        integer function adapt_X() result(ierr)
            use VMEC_vars, only: n_r_eq
            
            character(*), parameter :: rout_name = 'adapt_X'
            
            ! local variables
            character(len=max_str_ln) :: err_msg                                ! error message
            real(dp) :: one = 1.00000001_dp                                     ! one plus a little margin
            integer :: abs_min_n_X = 5                                          ! absolute minimum for the high-n theory
            integer :: rec_min_n_X = 10                                         ! recomended minimum for the high-n theory
            
            ! initialize ierr
            ierr = 0
            
            ! check n_X
            if (n_X.lt.abs_min_n_X) then
                n_X = abs_min_n_X
                call writo('WARNING: n_X has been increased to '//&
                    &trim(i2str(abs_min_n_X))//' but should be at least '//&
                    &trim(i2str(rec_min_n_X))//' for correct results')
            else if (n_X.lt.rec_min_n_X) then
                call writo('WARNING: n_X should be at least '//&
                    &trim(i2str(rec_min_n_X))//' for correct results')
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
        ! correct. max_m_X has to be greater  than min_m_X and nyq_fac has to be
        ! cat least 1
        integer function adapt_m() result (ierr)
            character(*), parameter :: rout_name = 'adapt_m'
            
            ! local variables
            character(len=max_str_ln) :: err_msg                                ! error message
            
            ! initialize ierr
            ierr = 0
            
            ! check min_m_X, max_m_X
            if (max_m_X.lt.min_m_X) then
                ierr = 1
                err_msg = 'max_m_X has to be larger or equal to min_m_X'
                CHCKERR(err_msg)
            end if
            
            ! adapt nyq_fac
            if (nyq_fac.lt.1) then
                call writo('WARNING: nyq_fac has been increased to 1')
                nyq_fac = 1
            end if
        end function adapt_m
        
        ! checks whether max_par - min_par indeed is a multiple of 2pi, which is
        ! a requisite for the correct functioning  of the code, as the integrals
        ! which are to be calculated in the parallel direction have to be closed
        ! loop intervals in the parallel direction
        integer function adapt_alpha() result(ierr)
            character(*), parameter :: rout_name = 'adapt_alpha'
            
            ! local variables
            real(dp) :: tol = 1E-5                                              ! tolerance
            real(dp) :: modulus
            character(len=max_str_ln) :: err_msg                                ! error message
            
            ! initialize ierr
            ierr = 0
            
            ! do the following checks only if calc_mesh_style = 0
            if (calc_mesh_style.eq.0) then
                ! check whether max_par > min_par
                if (max_par.lt.min_par+2*pi*(1-tol)) then
                    ierr = 1
                    err_msg = 'max_par has to be at least min_par + 2 pi'
                    CHCKERR(err_msg)
                end if
                
                ! calculate modulus
                modulus = mod(max_par-min_par,2*pi)
                modulus = min(modulus,2*pi-modulus)
                
                ! check whether tolerance is reached
                if (modulus/(2*pi).gt.tol) then
                    ierr = 1
                    err_msg = 'max_par - min_par has to be a multiple of 2 pi'
                    CHCKERR(err_msg)
                end if
            end if
            
            ! adapt n_alpha
            if (n_alpha.lt.1) then
                call writo('WARNING: n_alpha has been increased to 1')
                n_alpha = 1
            end if
        end function adapt_alpha
    end function read_input
end module input_ops
