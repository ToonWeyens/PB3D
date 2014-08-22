!------------------------------------------------------------------------------!
!   This module contains operations concerning giving input                    !
!------------------------------------------------------------------------------!
module input_ops
    use str_ops, only: strh2l
    use num_vars, only: &
        &dp, max_str_ln, style, min_alpha, max_alpha, n_alpha, max_it_NR, &
        &tol_NR, max_it_r, theta_var_along_B, input_i, n_seq_0, min_n_r_X, &
        &calc_mesh_style, EV_style, n_procs_per_alpha, plot_q, &
        &n_sol_requested, min_r_X, max_r_X
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
        &max_par, min_alpha, max_alpha, n_par, n_alpha, max_it_NR, &
        &tol_NR, max_it_r, theta_var_along_B, n_X, min_m_X, max_m_X, min_r_X, &
        &max_r_X, EV_style, n_procs_per_alpha, plot_q, n_sol_requested

contains
    ! queries for yes or no, depending on the flag yes:
    !   yes = .true.: yes is default answer
    !   yes = .false.: no is default answer
    logical function yes_no(yes)
        use output_ops, only: &
            &lvl_sep, lvl
        use time, only: start_time, stop_time
        
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
    subroutine read_input
        use num_vars, only: glob_rank
        
        ! local variables
        integer :: istat                                                        ! error
        
        if (glob_rank.eq.0) then                                                ! only global master
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
            calc_mesh_style = 0
            
            ! read user input
            if (input_i.ge.n_seq_0) then                                        ! otherwise, defaults are loaded
                read (input_i, nml=inputdata, iostat=istat)                     ! read input data
                if (istat.eq.0) then                                            ! input file succesfully read
                    call writo('Overwriting with user-provided file "' // &
                        &trim(input_name) // '"')
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
            
            max_it_NR = 50                                                      ! maximum 50 Newton-Rhapson iterations
            tol_NR = 1.0E-10_dp                                                 ! wanted relative error in Newton-Rhapson iteration
            max_it_r = 5                                                        ! maximum 5 levels of Richardson's extrapolation
            format_out = 1                                                      ! NETCDF output
            style = 1                                                           ! Richardson Extrapolation with normal discretization
            min_par = -4.0_dp*pi                                                ! minimum parallel angle
            max_par = 4.0_dp*pi                                                 ! maximum parallel angle
            n_par = 10                                                          ! number of parallel grid points
            min_alpha = 0.0_dp                                                  ! minimum field line label
            max_alpha = 2.0_dp*pi                                               ! maximum field line label
            n_alpha = 10                                                        ! number of different field lines
            theta_var_along_B = .true.                                          ! theta is used as the default parallel variable
            n_X = 20                                                            ! toroidal mode number of perturbation
            min_r_X = 0.1_dp                                                    ! minimum radius
            max_r_X = 1.0_dp                                                    ! maximum radius
            min_m_X = 20                                                        ! lowest poloidal mode number m_X
            max_m_X = 22                                                        ! highest poloidal mode number m_X
            EV_style = 1                                                        ! slepc solver for EV problem
            n_procs_per_alpha = 1                                               ! 1 processor per field line
            plot_q = .false.                                                    ! do not plot the q-profile with nq-m = 0
            n_sol_requested = 3                                                 ! request solutions with 3 highes EV
            min_n_r_X = 10                                                      ! at least 10 points in perturbation grid
        end subroutine
    end subroutine
end module input_ops
