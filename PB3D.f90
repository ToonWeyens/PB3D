! (from http://patorjk.com/software/taag/ ANSI shadow)

!   ██████╗ ██████╗ ██████╗ ██████╗ 
!   ██╔══██╗██╔══██╗╚════██╗██╔══██╗
!   ██████╔╝██████╔╝ █████╔╝██║  ██║
!   ██╔═══╝ ██╔══██╗ ╚═══██╗██║  ██║
!   ██║     ██████╔╝██████╔╝██████╔╝
!   ╚═╝     ╚═════╝ ╚═════╝ ╚═════╝ 

!------------------------------------------------------------------------------!
!>  Peeling Ballooning in 3D
!------------------------------------------------------------------------------!
!>  \author
!!  Toon Weyens,
!!  Contact: info@toonweyens.com
!------------------------------------------------------------------------------!
!>  \version    2.57
!!  \date       2012-2026
!!  Copyright (C) 2019-2026 Toon Weyens,
!!  Copyright (C) 2014-2019 ITER Organization
!!
!!  This program is free software: you can redistribute it and/or modify
!!  it under the terms of the GNU General Public License as published by
!!  the Free Software Foundation, either version 3 of the License, or
!!  (at your option) any later version.
!!
!!  This program is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!  GNU General Public License for more details.
!!
!!  You should have received a copy of the GNU General Public License
!!  along with this program. If not, see <https://www.gnu.org/licenses/>.
!------------------------------------------------------------------------------!
!>  \see
!!  References:
!!  \cite weyens2014theory
!!  \cite Weyens2017PB3D
!------------------------------------------------------------------------------!
#define CHCKERR if(ierr.ne.0) then; call sudden_stop(ierr); end if
program PB3D
    use num_vars, only: prog_name, prog_style, rank, rich_restart_lvl, eq_style
    use str_utilities, only: r2str, i2str
    use messages
    use grid_vars, only: grid_type
    use eq_vars, only: eq_1_type, eq_2_type
    use X_vars, only: X_1_type, X_2_type
    use vac_vars, only: vac_type
    use sol_vars, only: sol_type
    use HDF5_vars, only: init_HDF5
    use driver_eq, only: run_driver_eq
    use driver_X, only: run_driver_X
    use driver_sol, only: run_driver_sol
    use files_ops, only: open_input, open_output, parse_args, init_files, &
        &close_output
    use input_ops, only: read_input_opts, read_input_eq, print_output_in
    use input_utilities, only: dealloc_in
    use MPI_ops, only: start_MPI, stop_MPI, broadcast_input_opts, sudden_stop, &
        &stop_usage
    use MPI_utilities, only: broadcast_var
    use eq_ops, only: calc_normalization_const, normalize_input
    use eq_utilities, only: do_eq, eq_info
    use rich_ops, only: init_rich, term_rich, do_rich, start_rich_lvl, &
        &stop_rich_lvl
    use rich_vars, only: rich_info, &
        &rich_lvl
#if ldebug
    use num_vars, only: ltest
    use test, only: generic_tests
#endif
    
    implicit none

    ! local variables
    integer :: ierr                                                             ! error
    integer :: usage_status                                                     ! exit status after usage or help, -1 to go on
    type(grid_type), target :: grid_eq                                          ! equilibrium grid
    type(grid_type), pointer :: grid_eq_B => null()                             ! field-aligned equilibrium grid
    type(grid_type), target :: grid_X                                           ! perturbation grid
    type(grid_type), pointer :: grid_X_B => null()                              ! field-aligned perturbation grid
    type(grid_type) :: grid_sol                                                 ! solution grid
    type(eq_1_type) :: eq_1                                                     ! flux equilibrium variables
    type(eq_2_type) :: eq_2                                                     ! metric equilibrium variables
    type(X_1_type) :: X_1                                                       ! vectorial perturbation variables
    type(X_2_type) :: X_2                                                       ! integrated tensorial perturbation variables
    type(vac_type) :: vac                                                       ! vacuum variables
    type(sol_type) :: sol                                                       ! solution variables
    
    !------------------------------!
    !   Initialize some routines   !
    !------------------------------!
    prog_name = 'PB3D'                                                          ! program name
    ierr = start_MPI()                                                          ! start MPI
    CHCKERR
    prog_style = 1                                                              ! main part
    call print_hello()                                                          ! print message with time, etc
    call init_output()                                                          ! initialize output utilities
    call init_files()                                                           ! initialize file operations
    call init_time()                                                            ! initialize time
    call init_HDF5()                                                            ! initialize HDF5
 
    !-----------------------------------------------------------------------!
    !   Read the user-provided input file and the equilibrium code output   !
    !-----------------------------------------------------------------------!
    call start_time
    call writo('Initialization')
    call lvl_ud(1)
    usage_status = -1                                                           ! go on, unless rank 0 decides otherwise
    if (rank.eq.0) then                                                         ! only master
        ierr = parse_args(usage_status)                                         ! parse argument (options are used in open_input)
        CHCKERR
    end if
    ierr = broadcast_var(usage_status)                                          ! every process learns whether to stop
    CHCKERR
    if (usage_status.ge.0) call stop_usage(usage_status)                        ! all processes stop after usage or help
    if (rank.eq.0) then                                                         ! only master
        ierr = open_input()                                                     ! open the input files
        CHCKERR
        ierr = read_input_opts()                                                ! read input options file
        CHCKERR
        if (rich_restart_lvl.eq.1) then                                         ! start at first level
            ierr = read_input_eq()                                              ! read input equilibrium file
            CHCKERR
            ierr = calc_normalization_const()                                   ! set up normalization constants
            CHCKERR
            call normalize_input()                                              ! normalize the input
            ierr = open_output()                                                ! open output file
            CHCKERR
            ierr = print_output_in('in')                                        ! print input outputs
            CHCKERR
            call dealloc_in()                                                   ! clean up input from equilibrium codes
        else                                                                    ! restart at higher level
            ierr = open_output()                                                ! open output file
            CHCKERR
        end if
    end if
    ierr = broadcast_input_opts()                                               ! broadcast input options to other processors
    CHCKERR
    call writo('')
    call passed_time
    call writo('')
    call lvl_ud(-1)
    
#if ldebug
    !-------------------!
    !   Do some tests   !
    !-------------------!
    if (ltest) then
        call start_time
        call writo('Generic Tests')
        call lvl_ud(1)
        ierr = generic_tests()                                                  ! run generic test
        CHCKERR
        call writo('')
        call passed_time
        call writo('')
        call lvl_ud(-1)
    end if
#endif
    
    !----------------------------------------------!
    !   Initialize Richardson Extrapolation Loop   !
    !----------------------------------------------!
    call start_time
    ierr = init_rich()
    CHCKERR
    call stop_time
    
    RICH: do while(do_rich())
        !----------------------------!
        !   Start Richardson level   !
        !----------------------------!
        call start_time
        call writo('Richardson level '//trim(i2str(rich_lvl)))
        call lvl_ud(1)
        ierr = start_rich_lvl()                                                 ! start Richardson level, setting n_r_sol and other variables
        CHCKERR
        call stop_time
        call writo('')
        call lvl_ud(-1)
        
        PAR: do while(do_eq())
            !-----------------------------------!
            !   Main Driver: Equilibrium part   !
            !-----------------------------------!
            call start_time
            call writo('Equilibrium driver'//trim(rich_info())//&
                &trim(eq_info()))
            call lvl_ud(1)
            ierr = run_driver_eq(grid_eq,grid_eq_B,eq_1,eq_2,vac)               ! equilibrium driver
            CHCKERR
            call writo('')
            call passed_time
            call writo('')
            call lvl_ud(-1)
            
            !------------------------------------!
            !   Main driver: Perturbation part   !
            !------------------------------------!
            call start_time
            call writo('Perturbation driver'//trim(rich_info())//&
                &trim(eq_info()))
            call lvl_ud(1)
            ierr = run_driver_X(grid_eq,grid_eq_B,grid_X,grid_X_B,eq_1,eq_2,&
                &X_1,X_2)                                                       ! perturbation driver
            CHCKERR
            call writo('')
            call passed_time
            call writo('')
            call lvl_ud(-1)
        end do PAR
        
        ! For VMEC,  the metric equilibrium  variables are calculated  again in
        ! every equilibrium job, so the solution driver does not need those of
        ! the last job.
        if (eq_style.eq.1 .and. allocated(eq_2%jac_FD)) call eq_2%dealloc()
        
        !--------------------------------!
        !   Main driver: Solution part   !
        !--------------------------------!
        call start_time
        call writo('Solution driver'//trim(rich_info()))
        call lvl_ud(1)
        ierr = run_driver_sol(grid_eq,grid_X,grid_sol,X_2,vac,sol)              ! solution driver
        CHCKERR
        call writo('')
        call passed_time
        call writo('')
        call lvl_ud(-1)
        
        call start_time
        call stop_rich_lvl()
        call stop_time
    end do RICH
    
    !---------------------------------------!
    !   Stop Richarson Extrapolation Loop   !
    !---------------------------------------!
    call start_time
    call term_rich()                                                            ! stop Richardson loop
    call stop_time
    
    !--------------!
    !   Clean up   !
    !--------------!
    call writo('Clean up')
    call lvl_ud(1)
    if (eq_style.eq.2) call eq_2%dealloc()
    ierr = stop_MPI(grid_eq=grid_eq,&
        &grid_eq_B=grid_eq_B,&
        &grid_X=grid_X,&
        &grid_X_B=grid_X_B,&
        &grid_sol=grid_sol,&
        &eq_1=eq_1,&
        &eq_2=eq_2,&
        &X_1=X_1,&
        &X_2=X_2,&
        &vac=vac,&
        &sol=sol)
    CHCKERR
    call close_output
    call lvl_ud(-1)
    
    call print_goodbye
end program PB3D
