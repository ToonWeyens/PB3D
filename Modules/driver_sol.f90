!------------------------------------------------------------------------------!
!   Driver of the solution part of PB3D.                                       !
!------------------------------------------------------------------------------!
module driver_sol
#include <PB3D_macros.h>
    use str_ops
    use output_ops
    use messages
    use num_vars, only: dp, pi, max_str_ln
    use grid_vars, only: grid_type
    use eq_vars, only: eq_type
    use met_vars, only: met_type
    use X_vars, only: X_2_type
    use sol_vars, only: sol_type
    
    implicit none
    private
    public run_driver_sol
#if ldebug
    public debug_sol_grid
#endif
    
    ! global variables
#if ldebug
    logical :: debug_sol_grid = .false.                                         ! plot debug information for treatment of sol grid
    logical :: debug_X_norm = .true.                                           ! plot debug information X_norm
#endif
    
contains
    ! Implementation   of  the  driver,  using  Richardson's  extrapolation  and
    ! discretization in the normal direction
    ! [MPI] All ranks, parts are done  by global master only, other parts by the
    !       other ranks only
    !       Here, the tasks are split and allocated dynamically by the master to 
    !       the other groups
    integer function run_driver_sol() result(ierr)
        use num_vars, only: EV_style, rich_lvl_nr, max_it_r, n_sol_requested, &
            &no_guess, rank
        use grid_vars, only: dealloc_grid
        use eq_vars, only: dealloc_eq
        use met_vars, only: dealloc_met
        use X_vars, only: dealloc_X
        use sol_vars, only: dealloc_sol
        use utilities, only: test_max_memory
        use PB3D_ops, only: read_PB3D, reconstruct_PB3D
        use MPI_utilities, only: wait_MPI
        use SLEPC_ops, only: solve_EV_system_SLEPC
        use grid_ops, only: calc_norm_range, setup_and_calc_grid_sol
        use sol_ops, only: print_output_sol
        !!use utilities, only: calc_aux_utilities
#if ldebug
        use num_vars, only: iu, use_pol_flux_F
        use grid_ops, only: get_norm_interp_data
        use utilities, only: c, con
        use MPI_utilities, only: get_ser_var
        use eq_vars, only: vac_perm
#endif
        
        character(*), parameter :: rout_name = 'run_driver_sol'
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        type(grid_type) :: grid_eq, grid_sol                                    ! (field-aligned) equilibrium and solution grid
        type(eq_type) :: eq                                                     ! equilibrium variables
        type(met_type) :: met                                                   ! metric variables
        type(X_2_type) :: X                                                     ! tensorial perturbation variables
        type(sol_type) :: sol                                                   ! solution variables
        integer :: n_r_sol                                                      ! nr. of normal points in sol grid for Richardson lvl
        logical :: use_guess                                                    ! whether a guess is formed from previous level of Richardson
        logical :: done_richard                                                 ! is it converged?
        integer :: sol_limits(2)                                                ! min. and max. index of sol grid for this process
        real(dp), allocatable :: x_axis(:,:)                                    ! x axis for plot of Eigenvalues with n_r_sol
        real(dp), allocatable :: r_F_sol(:)                                     ! normal points in solution grid
        complex(dp), allocatable :: X_val_rich(:,:,:)                           ! Richardson array of eigenvalues X val
#if ldebug
        real(dp), allocatable :: loc_r_eq(:)                                    ! unrounded index of solution grid in equilibrium grid
        real(dp), allocatable :: ang_par_F(:,:,:)                               ! parallel angle theta_F or zeta_F, interpolated at sol grid
        real(dp), allocatable :: D2pkappa_nJ(:,:,:)                             ! D2p kappa_n J, interpolated at sol grid
        real(dp), allocatable :: Jh22mu_0(:,:,:)                                ! Jh22mu_0
        real(dp), allocatable :: q(:,:,:)                                       ! q
        complex(dp), allocatable :: PV_0(:,:,:)                                 ! PV_int_0, interpolated at sol grid
        complex(dp), allocatable :: X_norm(:,:,:)                               ! |X|^2 or other results to be plotted
        complex(dp), allocatable :: E_pot(:,:,:,:)                              ! E_pot to be plotted
        complex(dp), allocatable :: E_pot_int(:,:,:)                            ! B-averaged E_pot to be plotted
        integer :: i_lo, i_hi                                                   ! upper and lower index
        integer :: m, k                                                         ! counters
        integer :: id, jd, kd                                                   ! counters
        character(len=max_str_ln), allocatable :: var_names(:)                  ! names of variable to be plot
        character(len=max_str_ln) :: file_name                                  ! file name
#endif
        
        ! initialize ierr
        ierr = 0
        
        ! some preliminary things
        
        ! test maximum memory
        ierr = test_max_memory()
        CHCKERR('')
        
        !!! calculate auxiliary quantities for utilities
        !!call calc_aux_utilities                                                 ! calculate auxiliary quantities for utilities
        
        ! read PB3D output file
        ierr = read_PB3D(.false.,.true.,.false.,.true.,.false.)                 ! read the equilibrium and tensorial perturbation variables
        CHCKERR('')
        
        ! reconstruct PB3D variables
        ierr = reconstruct_PB3D(.false.,.true.,.false.,.true.,.false.,&
            &grid_eq=grid_eq,eq=eq,met=met,X_2=X)
        CHCKERR('')
        
        ! Initalize some variables for Richardson loop
        rich_lvl_nr = 1
        done_richard = .false.
        if (max_it_r.gt.1) then                                                 ! only do this if more than 1 Richardson level
            allocate(X_val_rich(1:max_it_r,1:max_it_r,1:n_sol_requested))
            X_val_rich = 0.0_dp
            allocate(x_axis(1:max_it_r,1:n_sol_requested))
        end if
        
        ! user output
        if (max_it_r.gt.1) then                                                 ! only do this if more than 1 Richardson level
            call writo('Starting perturbation calculation in field-aligned &
                &grid with Richardson extrapolation')
        else
            call writo('Starting perturbation calculation in field-aligned &
                &grid')
        end if
        call lvl_ud(1)                                                          ! before richardson loop
        
        Richard: do while (.not.done_richard .and. rich_lvl_nr.le.max_it_r)
            ! user output
            if (max_it_r.gt.1) then                                             ! only do this if more than 1 Richardson level
                call writo('Starting Level ' // trim(i2str(rich_lvl_nr)) // &
                    &' of Richardson extrapolation')
                call lvl_ud(1)                                                  ! beginning of one richardson loop
            end if
            
            ! Calculate number of  radial points for the  solution in Richardson
            ! loops and save in n_r_sol.
            call writo('calculating the normal points')
            call lvl_ud(1)
            call calc_n_r_sol(rich_lvl_nr,n_r_sol)
            call lvl_ud(-1)
            
            ! Divide solution grid under group processes, calculating the limits
            ! and the normal coordinate.
            allocate(r_F_sol(n_r_sol))
            ierr = calc_norm_range(sol_limits=sol_limits,r_F_sol=r_F_sol)
            CHCKERR('')
            
            ! create solution grid with division limits and setup normal
            ! coordinate
            call writo('Setting up solution grid')
            call lvl_ud(1)
            ierr = setup_and_calc_grid_sol(grid_eq,grid_sol,eq,r_F_sol,&
                &sol_limits)
            CHCKERR('')
            deallocate(r_F_sol)
            call lvl_ud(-1)
            
            ! set use_guess to .false. if no_guess
            if (no_guess) use_guess = .false.
            
            ! solve the system
            call writo('Solving the system')
            call lvl_ud(1)
            select case (EV_style)
                case(1)                                                         ! SLEPC solver for EV problem
                    ! solve the system
                    ierr = solve_EV_system_SLEPC(grid_eq,grid_sol,X,sol,&
                        &use_guess)
                    CHCKERR('')
                case default
                    err_msg = 'No EV solver style associated with '//&
                        &trim(i2str(EV_style))
                    ierr = 1
                    CHCKERR(err_msg)
            end select
            call lvl_ud(-1)
            call writo('System solved')
            
#if ldebug
            ! calculate |X|^2 directly from solution vector
            if (debug_X_norm) then
                call writo('Calculating |X|^2 directly from solution')
                call lvl_ud(1)
                
                ! allocate variables
                allocate(X_norm(grid_eq%n(1),grid_eq%n(2),grid_sol%loc_n_r))
                allocate(E_pot(grid_eq%n(1),grid_eq%n(2),grid_sol%loc_n_r,6))
                allocate(E_pot_int(grid_eq%n(2),grid_sol%loc_n_r,6))
                allocate(ang_par_F(grid_eq%n(1),grid_eq%n(2),grid_sol%loc_n_r))
                allocate(D2pkappa_nJ(grid_eq%n(1),grid_eq%n(2),grid_sol%loc_n_r))
                allocate(Jh22mu_0(grid_eq%n(1),grid_eq%n(2),grid_sol%loc_n_r))
                allocate(q(grid_eq%n(1),grid_eq%n(2),grid_sol%loc_n_r))
                allocate(PV_0(size(X%PV_int_0,1),grid_eq%n(2),grid_sol%loc_n_r))
                allocate(var_names(6))
                
                ! get normal interpolation data
                ierr = get_norm_interp_data(grid_eq,grid_sol,loc_r_eq)
                CHCKERR('')
                
                ! set variables
                do kd = 1,grid_sol%loc_n_r
                    i_lo = floor(loc_r_eq(kd))
                    i_hi = ceiling(loc_r_eq(kd))
                    if (use_pol_flux_F) then
                        ang_par_F(:,:,kd) = grid_eq%theta_F(:,:,i_lo) + &
                            &(grid_eq%theta_F(:,:,i_hi)-&
                            &grid_eq%theta_F(:,:,i_lo))*(loc_r_eq(kd)-i_lo)
                    else
                        ang_par_F(:,:,kd) = - (grid_eq%zeta_F(:,:,i_lo) + &
                            &(grid_eq%zeta_F(:,:,i_hi)-&
                            &grid_eq%zeta_F(:,:,i_lo))*(loc_r_eq(kd)-i_lo))
                    end if
                    !!!!!!!!!!!! TEMPORARILY !!!!!!!!!!
                    D2pkappa_nJ(:,:,kd) = &
                        &eq%pres_FD(i_lo,1)*eq%kappa_n(:,:,i_lo)*&
                        &met%jac_FD(:,:,i_lo,0,0,0) + &
                        &(eq%pres_FD(i_hi,1)*eq%kappa_n(:,:,i_hi)*&
                        &met%jac_FD(:,:,i_hi,0,0,0)-&
                        &eq%pres_FD(i_lo,1)*eq%kappa_n(:,:,i_lo)*&
                        &met%jac_FD(:,:,i_lo,0,0,0))*&
                        &(loc_r_eq(kd)-i_lo)
                    Jh22mu_0(:,:,kd) = (&
                        &met%h_FD(:,:,i_lo,4,0,0,0)*met%jac_FD(:,:,i_lo,0,0,0) + &
                        &(met%h_FD(:,:,i_hi,4,0,0,0)*met%jac_FD(:,:,i_hi,0,0,0)-&
                        &met%h_FD(:,:,i_lo,4,0,0,0)*met%jac_FD(:,:,i_lo,0,0,0))*&
                        &(loc_r_eq(kd)-i_lo))*vac_perm
                    q(:,:,kd) = &
                        &eq%q_saf_FD(i_lo,0) + (eq%q_saf_FD(i_hi,0)-&
                        &eq%q_saf_FD(i_lo,0))*(loc_r_eq(kd)-i_lo)
                    PV_0(:,:,kd) = X%PV_int_0(:,:,i_lo) + &
                        &(X%PV_int_0(:,:,i_hi)-&
                        &X%PV_int_0(:,:,i_lo))*(loc_r_eq(kd)-i_lo)
                    !!!!!!!!!!!! END TEMPORARILY !!!!!!!!!!
                end do
                X_norm = 0._dp
                
                ! calculate |X|^2
                var_names = 'X_norm'
                do kd = 1,grid_sol%loc_n_r
                    do m = 1,sol%n_mod
                        do k = 1,sol%n_mod
                            X_norm(:,:,kd) = X_norm(:,:,kd) + &
                                &conjg(sol%vec(k,kd,1))*sol%vec(m,kd,1)* &
                                &exp(iu*(X%m_1(k)-X%m_2(m))*ang_par_F(:,:,kd))
                        end do
                    end do
                end do
                
                ! plot |X|^2
                file_name = 'TEST_X_norm_PB3D'
                call plot_HDF5(var_names(1),file_name,realpart(X_norm),&
                    &tot_dim=[grid_eq%n(1:2),grid_sol%n(3)],&
                    &loc_offset=[0,0,grid_sol%i_min-1])
                
                !!!!!!!!!!!! TEMPORARILY !!!!!!!!!!
                ! Potential energies
                E_pot = 0._dp
                var_names = 'E_pot'
                do jd = 1,6
                    var_names(jd) = trim(var_names(jd))//'_'//trim(i2str(jd))
                end do
                
                ! calculate 1. 1/mu_0 J h22 |d/dtheta X|^2
                !!!!! ONLY FOR POL. FLUX !!!!
                do kd = 1,grid_sol%loc_n_r
                    do m = 1,sol%n_mod
                        do k = 1,sol%n_mod
                            E_pot(:,:,kd,1) = E_pot(:,:,kd,1) + &
                                &1/Jh22mu_0(:,:,kd)*&
                                &conjg(sol%vec(k,kd,1))*sol%vec(m,kd,1)*&
                                &(sol%n(m)*q(:,:,kd)-sol%m(m))*&
                                &(sol%n(k)*q(:,:,kd)-sol%m(k))*&
                                &exp(iu*(X%m_1(k)-X%m_2(m))*ang_par_F(:,:,kd))
                        end do
                    end do
                end do
                
                ! calculate 3. -2p' kappa_n |X|^2
                E_pot(:,:,:,3) = -2*D2pkappa_nJ*X_norm(:,:,:)
                
                ! plot E_pot
                file_name = 'TEST_E_pot_PB3D_RE'
                call plot_HDF5(var_names,file_name,realpart(E_pot),&
                    &tot_dim=[grid_eq%n(1:2),grid_sol%n(3),6],&
                    &loc_offset=[0,0,grid_sol%i_min-1,0])
                file_name = 'TEST_E_pot_PB3D_IM'
                call plot_HDF5(var_names,file_name,imagpart(E_pot),&
                    &tot_dim=[grid_eq%n(1:2),grid_sol%n(3),6],&
                    &loc_offset=[0,0,grid_sol%i_min-1,0])
                
                ! integrate along B (J is already built in!!)
                var_names = 'E_pot_int'
                E_pot_int = 0._dp
                do jd = 1,6
                    do id = 2,grid_eq%n(1)
                        E_pot_int(:,:,jd) = E_pot_int(:,:,jd) + &
                            &(E_pot(id,:,:,jd)+E_pot(id-1,:,:,jd))/2 &
                            &*(ang_par_F(id,:,:)-ang_par_F(id-1,:,:))
                    end do
                    var_names(jd) = trim(var_names(jd))//'_'//trim(i2str(jd))
                    
                    file_name = 'TEST_E_pot_PB3D_int_RE_'//trim(i2str(jd))
                    call print_GP_2D(var_names(jd),file_name,&
                        &realpart(transpose(E_pot_int(:,:,jd))),draw=.false.)
                    call draw_GP(var_names(jd),file_name,file_name,&
                        &grid_eq%n(2),1,.false.)
                    
                    file_name = 'TEST_E_pot_PB3D_int_IM_'//trim(i2str(jd))
                    call print_GP_2D(var_names(jd),file_name,&
                        &imagpart(transpose(E_pot_int(:,:,jd))),draw=.false.)
                    call draw_GP(var_names(jd),file_name,file_name,&
                        &grid_eq%n(2),1,.false.)
                end do
                
                var_names = 'E_pot_tot_int'
                file_name = 'TEST_E_pot_tot_PB3D_int_RE'
                call print_GP_2D(var_names(1),file_name,&
                    &realpart(transpose(sum(E_pot_int,3))),draw=.false.)
                call draw_GP(var_names(1),file_name,file_name,&
                    &grid_eq%n(2),1,.false.)
                
                file_name = 'TEST_E_pot_tot_PB3D_int_IM'
                call print_GP_2D(var_names(1),file_name,&
                    &imagpart(transpose(sum(E_pot_int,3))),draw=.false.)
                call draw_GP(var_names(1),file_name,file_name,&
                    &grid_eq%n(2),1,.false.)
                
                ! calculate X* PV_int_0 X (integrated along B) and save in 1
                E_pot_int = 0._dp
                !write(*,*) '    X = ', sol%vec(:,kd,1)
                do kd = 1,grid_sol%loc_n_r
                    !write(*,*) 'kd = ', kd
                    do m = 1,sol%n_mod
                        do k = 1,sol%n_mod
                            E_pot_int(:,kd,1) = E_pot_int(:,kd,1) + &
                                &conjg(sol%vec(k,kd,1))*sol%vec(m,kd,1)* &
                                &con(PV_0(c([k,m],.true.,n=X%n_mod(1)),&
                                &:,kd),[k,m],.true.,[grid_eq%n(2)])
                        !!write(*,*) 'kd,k,m',kd,k,m, 'using', &
                                !!&con(PV_0(c([k,m],.true.,n=X%n_mod(1)),&
                                !!&:,kd),[k,m],.true.,[grid_eq%n(2)]), &
                                !!&' result',  &
                                !!&conjg(sol%vec(k,kd,1))*sol%vec(m,kd,1)* &
                                !!&con(PV_0(c([k,m],.true.,n=X%n_mod(1)),&
                                !!&:,kd),[k,m],.true.,[grid_eq%n(2)])
                        !!read(*,*)
                        !write(*,*) '        PV(',k,m,') =', &
                                !&con(PV_0(c([k,m],.true.,n=X%n_mod(1)),&
                                !&:,kd),[k,m],.true.,[grid_eq%n(2)])
                        end do
                    end do
                    !write(*,*) '    E_pot_int = '//trim(c2str(E_pot_int(1,kd,1)))
                end do
                
                file_name = 'TEST_E_pot_tot_MAT_PB3D_int_RE'
                call print_GP_2D(var_names(1),file_name,&
                    &realpart(transpose(E_pot_int(:,:,1))),draw=.false.)
                call draw_GP(var_names(1),file_name,file_name,&
                    &grid_eq%n(2),1,.false.)
                
                file_name = 'TEST_E_pot_tot_MAT_PB3D_int_IM'
                call print_GP_2D(var_names(1),file_name,&
                    &imagpart(transpose(E_pot_int(:,:,1))),draw=.false.)
                call draw_GP(var_names(1),file_name,file_name,&
                    &grid_eq%n(2),1,.false.)
                
                ! integrate in normal direction
                call writo('Crude integral: '//&
                    &trim(c2str(sum(E_pot_int(:,:,1))*&
                    &(grid_sol%r_F(grid_sol%n(3))-grid_sol%r_F(1))/&
                    &(grid_sol%n(3)-1)))//' with step size '//r2str(&
                    &(grid_sol%r_F(grid_sol%n(3))-grid_sol%r_F(1))/&
                    &(grid_sol%n(3)-1)))
                !!!!!!!!!!!! END TEMPORARILY !!!!!!!!!!
                
                call writo('The output should be compared with the POST-&
                    &output')
                
                call lvl_ud(-1)
            end if
#endif
            
            ! write solution variables to output
            ierr = print_output_sol(grid_sol,sol)
            CHCKERR('')
            
            ! Richardson extrapolation
            if (max_it_r.gt.1) then                                             ! only do this if more than 1 Richardson level
                ! update the x axis of the Eigenvalue plot
                x_axis(rich_lvl_nr,:) = 1.0_dp*n_r_sol
                
                ! update  the  variable  X_val_rich  with the  results  of  this
                ! Richardson level
                call writo('updating Richardson extrapolation variables')
                call lvl_ud(1)
                ierr = calc_rich_ex(rich_lvl_nr,sol%val,X_val_rich,&
                    &done_richard,use_guess)
                CHCKERR('')
                call lvl_ud(-1)
            else                                                                ! if not, Richardson is done
                done_richard = .true.
            end if
            
            ! user output
            if (max_it_r.gt.1) then                                             ! only do this if more than 1 Richardson level
                call lvl_ud(-1)                                                 ! beginning of one richardson loop
                call writo('Completed Level ' // trim(i2str(rich_lvl_nr)) // &
                    &' of Richardson extrapolation')
            end if
            
            ! draw output
            if (max_it_r.gt.1 .and. rank.eq.0) call draw_X_val_rich(X_val_rich)
            
            ! clean up
            call dealloc_eq(eq)
            call dealloc_met(met)
            call dealloc_grid(grid_eq)
            call dealloc_grid(grid_sol)
            
            ! synchronize MPI
            ierr = wait_MPI()
            CHCKERR('')
        end do Richard
        
        ! user output
        call lvl_ud(-1)                                                         ! done with richardson
        if (max_it_r.gt.1) then                                                 ! only do this if more than 1 Richardson level
            call writo('Finished Richardson loop')
        else
            call writo('Finished perturbation calculation')
        end if
        
        ! clean up
        call dealloc_X(X)
        call dealloc_sol(sol)
    contains
        ! Calculates the number  of normal points n_r_sol for  the solution grid
        ! for the various Richardson iterations.
        ! The aim is to  halve the step size, which is given  by dx(n) = 1/(n-1)
        ! or, inverting: n(dx) = 1 + 1/dx.
        ! This yields n(dx/2)/n(dx) = (2+dx)/(1+dx) = (2n(dx)-1)/n(dx)
        ! The recursion formula is therefore: n(dx/2) = 2n(dx) - 1
        subroutine calc_n_r_sol(ir,n_r_sol)
            use X_vars, only: min_n_r_sol
            
            ! input / output
            integer, intent(in) :: ir
            integer, intent(inout) :: n_r_sol
            
            write(*,*) '!!!! TEMPORARILY DISABLED BECAUSE RICH. EXT. NOT &
                &WORKING PROPERLY !!!'
            
            n_r_sol = min_n_r_sol*ir
            !!!if (ir.eq.1) then
                !!!n_r_sol = min_n_r_sol
            !!!else
                !!!n_r_sol = 2 * n_r_sol - 1
            !!!end if
            call writo(trim(i2str(n_r_sol))//' normal points for this level')
        end subroutine calc_n_r_sol
        
        ! Calculates  the  coefficients of  the  Eigenvalues  in the  Richardson
        ! extrapolation.
        ! This is done using the recursive formula
        !   X_val_rich(ir,ir2,:) = X_val_rich(ir,ir2-1,:) +  1/(2^(2ir2) - 1) * 
        !       (X_val_rich(ir,ir2-1,:) - X_val_rich(ir-1,ir2-1,:)),
        ! as described in [ADD REF]
        ! [MPI] All ranks
        integer function calc_rich_ex(ir,X_val,X_val_rich,done_richard,&
                &use_guess_for_next_level) result(ierr)
            use num_vars, only: tol_r
            
            character(*), parameter :: rout_name = 'calc_rich_ex'
            
            ! input / output
            integer, intent(inout) :: ir                                        ! level of Richardson extrapolation (starting at 1)
            complex(dp), intent(in) :: X_val(:)                                 ! EV for this Richardson level
            complex(dp), intent(inout) :: X_val_rich(:,:,:)                     ! extrapolated coefficients
            logical, intent(inout) :: done_richard                              ! if Richardson loop has converged sufficiently
            logical, intent(inout) :: use_guess_for_next_level                  ! if a guessed is used for next Richardson level
            
            ! local variables
            character(len=max_str_ln) :: err_msg                                ! error message
            integer :: ir2                                                      ! counter
            complex(dp), allocatable :: corr(:)                                 ! correction
            real(dp) :: max_corr                                                ! maximum of maximum of correction
            real(dp), save :: prev_max_corr                                     ! max_corr at previous Richardson level
            integer :: loc_max_corr(1)                                          ! location of maximum of correction
            
            ! initialize ierr
            ierr = 0
            write(*,*) '!!! FOR RICHARDSON EXTRAPOLATION YOU HAVE TO FIND &
                &CONSTANT EIGENVALUES !!!'
            
            ! tests
            if (size(X_val_rich,1).ne.size(X_val_rich,2) .or. &
                &ir.gt.size(X_val_rich,1)) then
                ierr = 1
                err_msg = 'X_val_rich has to have correct dimensions'
                CHCKERR(err_msg)
            end if
            
            ! allocate correction
            allocate(corr(size(X_val))); corr = 1.0E15
            
            ! do calculations if ir > 1
            X_val_rich(ir,1,1:size(X_val)) = X_val
            do ir2 = 2,ir
                corr = 1./(2**(2*ir2)-1.) * &
                    &(X_val_rich(ir,ir2-1,1:size(X_val)) - &
                    &X_val_rich(ir-1,ir2-1,1:size(X_val)))
                X_val_rich(ir,ir2,1:size(X_val)) = &
                    &X_val_rich(ir,ir2-1,1:size(X_val)) + corr
            end do
            
            if (ir.gt.1) then                                                   ! only do this if in Richardson level higher than 1
                ! get maximum and location of maximum for relative correction
                max_corr = maxval(abs(corr/X_val_rich(ir,ir,1:size(X_val))))
                loc_max_corr = maxloc(abs(corr/X_val_rich(ir,ir,1:size(X_val))))
                call writo('maximum relative error: '//trim(r2strt(max_corr))//&
                    &' for Eigenvalue '//trim(i2str(loc_max_corr(1))))
                
                ! check whether tolerance has been  reached for last value ir2 =
                ! ir
                if (maxval(abs(corr/X_val_rich(ir,ir,1:size(X_val))))&
                    &.lt.tol_r) then
                    done_richard = .true.
                    call writo('tolerance '//trim(r2strt(tol_r))//&
                        &' reached after '//trim(i2str(ir))//' iterations')
                else
                    call writo('tolerance '//trim(r2strt(tol_r))//' not yet &
                        &reached')
                end if
                
                ! determine use_guess_for_next_level
                use_guess_for_next_level = max_corr.lt.prev_max_corr            ! set guess for next level to true if max_corr is decreasing
                prev_max_corr = max_corr                                        ! save max_corr for next level
            else
                use_guess_for_next_level = .true.                               ! for first Richardson level, set guess for next level to true
            end if
            
            !!! TEMPORARILY !!!
            call writo('!!!!! DO NOT USE GUESS FOR MODIFIED RICHARDSON !!!')
            use_guess_for_next_level = .false.
            
            ! check for convergence
            if (.not.done_richard) then
                if (ir.lt.max_it_r) then                                        ! not yet at maximum Richardson iteration
                    ir = ir + 1
                else                                                            ! maximum nr. of Richardson iterations reached
                    call writo('maximum number of Richardson iterations &
                        &reached')
                    done_richard = .true.
                end if
            end if
        end function calc_rich_ex
        
        ! Draws  the  Eigenvalues   for   the  different  levels  of  Richardson
        ! extrapolation as a function of the number of normal points.
        subroutine draw_X_val_rich(X_val_rich)
            ! input / output
            complex(dp), allocatable :: X_val_rich(:,:,:)                       ! Richardson array of eigenvalues X val
            
            ! local variables
            character(len=max_str_ln) :: plot_title                             ! title for plots
            character(len=max_str_ln) :: plot_name                              ! name of plot
            character(len=max_str_ln) :: draw_ops                               ! optional drawing options
            
            ! user output
            call writo('Plotting Eigenvalues as function of nr. of normal &
                &points in Richardson Extrapolation')
            call lvl_ud(1)
            
            ! set up drawing options
            draw_ops = 'pt 7 ps 0.2'
            
            ! output on screen
            plot_title = 'Eigenvalues as function of nr. of normal points'
            plot_name = 'Eigenvalues_richardson'
            call print_GP_2D(plot_title,plot_name,&
                &realpart(X_val_rich(1:rich_lvl_nr,1,:)),&
                &x=x_axis(1:rich_lvl_nr,:),draw=.false.)
            
            ! same output in file as well
            call draw_GP(plot_title,plot_name,plot_name,&
                &n_sol_requested,1,.false.,draw_ops=draw_ops)
            
            ! user output
            call lvl_ud(-1)
            call writo('Done plotting Eigenvalues')
        end subroutine draw_X_val_rich
    end function run_driver_sol
end module driver_sol
