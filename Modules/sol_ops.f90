!------------------------------------------------------------------------------!
!   Operations  on the  solution  vectors  such as  returning  X,  U or  their !
!   derivates, etc.                                                            !
!------------------------------------------------------------------------------!
module sol_ops
#include <PB3D_macros.h>
    use str_ops
    use output_ops
    use messages
    use num_vars, only: dp, iu, max_str_ln, pi
    use grid_vars, only: grid_type
    use eq_vars, only: eq_type
    use met_vars, only: met_type
    use X_vars, only: X_1_type
    use sol_vars, only: sol_type

    implicit none
    private
    public calc_XUQ, plot_X_vec, decompose_energy, print_output_sol, plot_X_vals
#if ldebug
    public debug_calc_XUQ_arr, debug_plot_X_vec, debug_calc_E, debug_DU
#endif
    
    ! interfaces
    interface calc_XUQ
        module procedure calc_XUQ_arr, calc_XUQ_ind
    end interface
    
    ! global variables
#if ldebug
    logical :: debug_calc_XUQ_arr = .false.                                     ! plot debug information for calc_XUQ_arr
    logical :: debug_plot_X_vec = .false.                                       ! plot debug information for plot_X_vec
    logical :: debug_calc_E = .false.                                           ! plot debug information for calc_E
    logical :: debug_X_norm = .false.                                           ! plot debug information X_norm
    logical :: debug_DU = .false.                                               ! plot debug information for calculation of DU
#endif

contains
    ! Calculates the normal or geodesic components of the plasma perturbation or
    ! the magnetic perturbation, indicated by the variable "XUQ_style":
    !   - XUQ_style = 1: X (supports parallel derivative)
    !   - XUQ_style = 2: U (supports parallel derivative)
    !   - XUQ_style = 3: Qn
    !   - XUQ_style = 4: Qg
    ! X and U support the calculation  of the parallel derivative and for Qn and
    ! Qg, the metric variables have to be provided as well.
    ! As  this   operation  requires   the  equilibrium,   metric,  perturbation
    ! and  solution  variables,  the  grid  in which  this  happens  is  not  so
    ! straightforward:
    ! The equilibrium and metric variables are tabulated in an equilibrium grid,
    ! which needs to  have the same angular extent as  the perturbation grid, in
    ! which the perturbation variables are tabulated. The normal extent of these
    ! two grids can  be different, though, as  well as the normal  extent of the
    ! solution grid, in which the solution variables are tabulated. 
    ! Currently,  the  perturbation variables  are  supposed  to have  the  same
    ! angular  extent as  the equilibrium  and  metric variables,  and the  same
    ! normal  extent as  the solution  variables.  The output  then follows  the
    ! perturbation grid.
    integer function calc_XUQ_arr(grid_eq,grid_X,grid_sol,eq,X,sol,X_id,&
        &XUQ_style,time,XUQ,met,deriv) result(ierr)                             ! (time) array version
        use num_vars, only: use_pol_flux_F, norm_disc_prec_sol
        use utilities, only: con2dis, calc_deriv
        use grid_ops, only: get_norm_interp_data
#if ldebug
        use utilities, only: calc_int
#endif
        
        character(*), parameter :: rout_name = 'calc_XUQ_arr'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid
        type(grid_type), intent(in) :: grid_X                                   ! perturbation grid
        type(grid_type), intent(in) :: grid_sol                                 ! solution grid
        type(eq_type), intent(in) :: eq                                         ! equilibrium variables
        type(X_1_type), intent(in) :: X                                         ! perturbation variables
        type(sol_type), intent(in) :: sol                                       ! solution variables
        integer, intent(in) :: X_id                                             ! nr. of Eigenvalue
        integer, intent(in) :: XUQ_style                                        ! whether to calculate X, U, Qn or Qg
        real(dp), intent(in) :: time(:)                                         ! time range
        complex(dp), intent(inout) :: XUQ(:,:,:,:)                              ! X, U, Qn or Qg
        type(met_type), optional, intent(in) :: met                             ! metric variables
        logical, intent(in), optional :: deriv                                  ! return parallel derivative
        
        ! local variables
        integer :: n_t                                                          ! number of time points
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer :: id, jd, kd                                                   ! counter
        complex(dp) :: sqrt_X_val_norm                                          ! normalized sqrt(X_val)
        logical :: deriv_loc                                                    ! local copy of deriv
        complex(dp), allocatable :: DX_vec(:)                                   ! normal derivative of X_vec for a specific mode
        complex(dp), allocatable :: fac_0(:,:,:), fac_1(:,:,:)                  ! factor to be multiplied with X and DX
        complex(dp), allocatable :: XUQ_loc(:,:)                                ! complex XUQ without time at a normal point
        real(dp), allocatable :: par_fac(:)                                     ! multiplicative factor due to parallel derivative, without iu
        real(dp), allocatable :: expon(:,:,:)                                   ! exponent
        real(dp), allocatable :: loc_r_eq(:)                                    ! loc_r_F of sol grid interpolated in eq grid
        real(dp), allocatable :: jq(:)                                          ! iota or q, interpolated at solution grid
        real(dp), allocatable :: S(:,:,:), inv_J(:,:,:)                         ! S and 1/J, interpolated at solution grid
        integer :: i_lo, i_hi                                                   ! upper and lower index for interpolation of eq grid to sol grid
#if ldebug
        real(dp), allocatable :: X_vec_ALT(:)                                   ! X_vec calculated from DX_vec
#endif
        
        ! initialize ierr
        ierr = 0
        
        ! set up n_t
        n_t = size(time)
        
        ! tests
        if (grid_eq%n(1).ne.grid_X%n(1) .or. grid_eq%n(2).ne.grid_X%n(2) .or. &
            &grid_sol%loc_n_r.ne.grid_X%loc_n_r) then
            ierr = 1
            err_msg = 'Grids need to be compatible'
            CHCKERR(err_msg)
        end if
        if (grid_eq%n(1).ne.size(XUQ,1) .or. grid_eq%n(2).ne.size(XUQ,2) .or. &
            &grid_X%loc_n_r.ne.size(XUQ,3) .or. n_t.ne.size(XUQ,4)) then
            ierr = 1
            err_msg = 'XUQ needs to have the correct dimensions'
            CHCKERR(err_msg)
        end if
        if ((XUQ_style.eq.3 .or. XUQ_style.eq.4).and..not.present(met)) then
            ierr = 1
            err_msg = 'When calculating Q, also metric variables have to be &
                &provided'
            CHCKERR(err_msg)
        end if
        
        ! set up local copy of deriv
        deriv_loc = .false.
        if (present(deriv)) deriv_loc = deriv
        
        ! set up normalized sqrt(X_val)
        sqrt_X_val_norm = sqrt(sol%val(X_id))
        if (imagpart(sqrt_X_val_norm).gt.0) &
            &sqrt_X_val_norm = - sqrt_X_val_norm                                ! exploding solution, not the decaying one
        sqrt_X_val_norm = sqrt_X_val_norm / abs(sqrt_X_val_norm)                ! normalize
        
        ! allocate helper variables
        ! Multiplicative factors fac_0 and fac_1 and par_fac are used in
        ! XUQ = fac_0*X + fac_1*DX with fac_0 and fac_1:
        !   1. fac_0 = 1 (i(nq-m) or i(n-iota m) for deriv)
        !      fac_1 = 0
        !   2. fac_0 = U_0 (DU_0 for deriv)
        !      fac_1 = U_1 (DU_1 for deriv)
        !   3. fac_0 = i(nq-m)/J or i(n-iota m)/J (no deriv)
        !      fac_1 = 0
        !   4. fac_0 = DU_0/J - S (no deriv)
        !      fac_1 = DU_1/J
        ! where par_fac = nq-m or n-iota m is used.
        ! Note that if the resolution for X_vec is bad, the numerical derivative
        ! is very  inaccurate, therefore  only smooth (i.e.  physical) solutions
        ! should be looked at.
        allocate(DX_vec(grid_X%loc_n_r))
        allocate(XUQ_loc(size(XUQ,1),size(XUQ,2)))
        allocate(fac_0(grid_X%n(1),grid_X%n(2),grid_X%loc_n_r))
        allocate(fac_1(grid_X%n(1),grid_X%n(2),grid_X%loc_n_r))
        allocate(par_fac(grid_X%loc_n_r))
        allocate(jq(grid_X%loc_n_r))
        allocate(S(grid_X%n(1),grid_X%n(2),grid_X%loc_n_r))
        if (present(met)) allocate(inv_J(grid_X%n(1),grid_X%n(2),&
            &grid_X%loc_n_r))
        allocate(expon(grid_X%n(1),grid_X%n(2),grid_X%loc_n_r))
        
        ! calculate  normal  interpolation  tables for  equilibrium  grid  (this
        ! concerns eq variables and met variables)
        ierr = get_norm_interp_data(grid_eq,grid_X,loc_r_eq)
        CHCKERR('')
        
        ! interpolate iota or q, S and, if present, 1/J
        do kd = 1,grid_X%loc_n_r
            ! set lower and upper index for this normal point
            i_lo = floor(loc_r_eq(kd))
            i_hi = ceiling(loc_r_eq(kd))
            
            if (use_pol_flux_F) then
                jq(kd) = eq%q_saf_FD(i_lo,0) + (loc_r_eq(kd)-i_lo) * &
                    &(eq%q_saf_FD(i_hi,0)-eq%q_saf_FD(i_lo,0))
            else
                jq(kd) = eq%rot_t_FD(i_lo,0) + (loc_r_eq(kd)-i_lo) * &
                    &(eq%rot_t_FD(i_hi,0)-eq%rot_t_FD(i_lo,0))
            end if
            S(:,:,kd) = eq%S(:,:,i_lo) + (loc_r_eq(kd)-i_lo) * &
                &(eq%S(:,:,i_hi)-eq%S(:,:,i_lo))
            if (present(met)) inv_J(:,:,kd) = &
                &1._dp/met%jac_FD(:,:,i_lo,0,0,0) + (loc_r_eq(kd)-i_lo) * &
                &(1._dp/met%jac_FD(:,:,i_hi,0,0,0)-&
                &1._dp/met%jac_FD(:,:,i_lo,0,0,0))
        end do
        
        ! initialize XUQ
        XUQ = 0._dp
        
        ! iterate over all Fourier modes
        Fourier: do jd = 1,sol%n_mod
            ! set up parallel multiplicative factor par_fac in normal eq grid
            if (use_pol_flux_F) then
                par_fac = sol%n(jd)*jq-sol%m(jd)
            else
                par_fac = sol%n(jd)-sol%m(jd)*jq
            end if
            
            ! set up exponent
            expon = sol%n(jd)*grid_X%zeta_F - sol%m(jd)*grid_X%theta_F
            
            ! set up multiplicative factors in eq grid
            select case (XUQ_style)
                case (1)                                                        ! calculating X
                    if (deriv_loc) then                                         ! parallel derivative
                        do kd = 1,grid_X%loc_n_r
                            fac_0(:,:,kd) = iu*par_fac(kd)
                        end do
                    else
                        fac_0 = 1._dp
                    end if
                    fac_1 = 0._dp
                case (2)                                                        ! calculating U
                    if (deriv_loc) then                                         ! parallel derivative
                        fac_0 = X%DU_0(:,:,:,jd)
                        fac_1 = X%DU_1(:,:,:,jd)
                    else
                        fac_0 = X%U_0(:,:,:,jd)
                        fac_1 = X%U_1(:,:,:,jd)
                    end if
                case (3)                                                        ! calculating Qn
                    if (deriv_loc) then                                         ! parallel derivative
                        ierr = 1
                        err_msg = 'For XUQ_style = '//trim(i2str(XUQ_style))//&
                            &' calculation of parallel derivative not supported'
                        CHCKERR('')
                    else
                        do kd = 1,grid_X%loc_n_r
                            fac_0(:,:,kd) = iu*par_fac(kd)*inv_J(:,:,kd)
                        end do
                    end if
                    fac_1 = 0._dp
                case (4)                                                        ! calculating Qg
                    if (deriv_loc) then                                         ! parallel derivative
                        ierr = 1
                        err_msg = 'For XUQ_style = '//trim(i2str(XUQ_style))//&
                            &' calculation of parallel derivative not supported'
                        CHCKERR('')
                    else
                        do kd = 1,grid_X%loc_n_r
                            fac_0(:,:,kd) = X%DU_0(:,:,kd,jd)*inv_J(:,:,kd) &
                                &- S(:,:,kd)
                            fac_1(:,:,kd) = X%DU_1(:,:,kd,jd)*inv_J(:,:,kd)
                        end do
                    end if
                case default
                    err_msg = 'No XUQ style associated with '//&
                        &trim(i2str(XUQ_style))
                    ierr = 1
                    CHCKERR(err_msg)
            end select
            
            ! set up normal derivative of X vec
            if (XUQ_style.eq.2 .or. XUQ_style.eq.4) then
                ierr = calc_deriv(sol%vec(jd,:,X_id),DX_vec,grid_X%loc_r_F,1,&
                    &norm_disc_prec_sol)
                CHCKERR('')
#if ldebug
                if (debug_calc_XUQ_arr) then
                    call writo('For mode '//trim(i2str(jd))//', testing &
                        &whether DX_vec is correct by comparing its integral &
                        &with original X_vec')
                    allocate(X_vec_ALT(1:grid_X%loc_n_r))
                    ierr = calc_int(realpart(DX_vec),&
                        &grid_X%loc_r_F(1:grid_X%loc_n_r),X_vec_ALT)
                    CHCKERR('')
                    call print_GP_2D('TEST_RE_DX_vec','',reshape(&
                        &[realpart(sol%vec(jd,:,X_id)),X_vec_ALT],&
                        &[grid_X%loc_n_r,2]),x=reshape(&
                        &[grid_X%loc_r_F(1:grid_X%loc_n_r),&
                        &grid_X%loc_r_F(1:grid_X%loc_n_r)],&
                        &[grid_X%loc_n_r,2]))
                    ierr = calc_int(imagpart(DX_vec),&
                        &grid_X%loc_r_F(1:grid_X%loc_n_r),X_vec_ALT)
                    CHCKERR('')
                    call print_GP_2D('TEST_IM_DX_vec','',reshape(&
                        &[imagpart(sol%vec(jd,:,X_id)),X_vec_ALT],&
                        &[grid_X%loc_n_r,2]),x=reshape(&
                        &[grid_X%loc_r_F(1:grid_X%loc_n_r),&
                        &grid_X%loc_r_F(1:grid_X%loc_n_r)],&
                        &[grid_X%loc_n_r,2]))
                    deallocate(X_vec_ALT)
                end if
#endif
            else
                DX_vec = 0._dp
            end if
            
            ! iterate over all normal points in sol grid
            normal: do kd = 1,grid_X%loc_n_r
                ! set lower and upper index for this normal point
                i_lo = floor(loc_r_eq(kd))
                i_hi = ceiling(loc_r_eq(kd))
                
                ! set up loc complex XUQ without time at this normal point
                XUQ_loc = exp(iu*expon(:,:,kd)) * &
                    &(sol%vec(jd,kd,X_id)*fac_0(:,:,kd) + &
                    &DX_vec(kd)*fac_1(:,:,kd))
                
                ! iterate over time steps
                do id = 1,n_t
                    ! add current mode
                    XUQ(:,:,kd,id) = XUQ(:,:,kd,id) + XUQ_loc * &
                        &exp(iu*sqrt_X_val_norm*time(id)*2*pi)
                end do
            end do normal
        end do Fourier
    end function calc_XUQ_arr
    integer function calc_XUQ_ind(grid_eq,grid_X,grid_sol,eq,X,sol,X_id,&
        &XUQ_style,time,XUQ,met,deriv) result(ierr)                             ! (time) individual version
        character(*), parameter :: rout_name = 'calc_XUQ_ind'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! equilibirum grid
        type(grid_type), intent(in) :: grid_X                                   ! perturbation grid
        type(grid_type), intent(in) :: grid_sol                                 ! solution grid
        type(eq_type), intent(in) :: eq                                         ! equilibrium variables
        type(X_1_type), intent(in) :: X                                         ! perturbation variables
        type(sol_type), intent(in) :: sol                                       ! solution variables
        integer, intent(in) :: X_id                                             ! nr. of Eigenvalue
        integer, intent(in) :: XUQ_style                                        ! whether to calculate X, U, Qn or Qg
        real(dp), intent(in) :: time                                            ! time
        complex(dp), intent(inout) :: XUQ(:,:,:)                                ! normal component of perturbation
        type(met_type), optional, intent(in) :: met                             ! metric variables
        logical, intent(in), optional :: deriv                                  ! return parallel derivative
        
        ! local variables
        complex(dp), allocatable :: XUQ_arr(:,:,:,:)                            ! dummy variable to use array version
        
        ! allocate array version of XUQ
        allocate(XUQ_arr(size(XUQ,1),size(XUQ,2),size(XUQ,3),1))
        
        ! call array version
        ierr = calc_XUQ_arr(grid_eq,grid_X,grid_sol,eq,X,sol,X_id,XUQ_style,&
            &[time],XUQ_arr,met,deriv)
        CHCKERR('')
        
        ! copy array to individual XUQ
        XUQ = XUQ_arr(:,:,:,1)
        
        ! deallocate array version of XUQ
        deallocate(XUQ_arr)
    end function calc_XUQ_ind
    
    ! plots Eigenvalues
    integer function plot_X_vals(sol,last_unstable_id) result(ierr)
        use num_vars, only: rank
        use MPI_utilities, only: wait_MPI
        
        character(*), parameter :: rout_name = 'plot_X_vals'
        
        ! input / output
        type(sol_type), intent(in) :: sol                                       ! solution variables
        integer, intent(in) :: last_unstable_id                                 ! index of last unstable EV
        
        ! local variables
        character(len=max_str_ln) :: plot_title                                 ! title for plots
        character(len=max_str_ln) :: plot_name                                  ! file name for plots
        integer :: n_sol_found                                                  ! how many solutions found and saved
        integer :: id                                                           ! counter
        
        ! initialize ierr
        ierr = 0
        
        ! set local variables
        n_sol_found = size(sol%val)
        
        ! only let master plot
        if (rank.eq.0) then
            ! Last Eigenvalues
            plot_title = 'final Eigenvalues omega^2 [log]'
            plot_name = 'Eigenvalues'
            call print_GP_2D(plot_title,plot_name,&
                &log10(abs(realpart(sol%val(1:n_sol_found)))),draw=.false.)
            call draw_GP(plot_title,plot_name,plot_name,1,1,.false.)
            
            ! Last Eigenvalues: unstable range
            if (last_unstable_id.gt.0) then
                plot_title = 'final unstable Eigenvalues omega^2'
                plot_name = 'Eigenvalues_unstable'
                call print_GP_2D(plot_title,plot_name,&
                    &realpart(sol%val(1:last_unstable_id)),&
                    &x=[(id*1._dp,id=1,last_unstable_id)],draw=.false.)
                call draw_GP(plot_title,plot_name,plot_name,1,1,.false.)
            end if
            
            ! Last Eigenvalues: stable range
            if (last_unstable_id.lt.n_sol_found) then
                plot_title = 'final stable Eigenvalues omega^2'
                plot_name = 'Eigenvalues_stable'
                call print_GP_2D(plot_title,plot_name,&
                    &realpart(sol%val(last_unstable_id+1:n_sol_found)),&
                    &x=[(id*1._dp,id=last_unstable_id+1,n_sol_found)],&
                    &draw=.false.)
                call draw_GP(plot_title,plot_name,plot_name,1,1,.false.)
            end if
        end if
        
        ! synchronize processes
        ierr = wait_MPI()
        CHCKERR('')
    end function plot_X_vals
    
    ! Plots  Eigenvectors  using  the  angular  part  of  the  the  provided
    ! equilibrium  grid and  the normal  part of  the provided  perturbation
    ! grid.
#if ldebug
    integer function plot_X_vec(grid_eq,grid_X,grid_sol,eq,met,X,sol,XYZ,X_id,&
        &res_surf) result(ierr)
#else
    integer function plot_X_vec(grid_eq,grid_X,grid_sol,eq,X,sol,XYZ,X_id,&
        &res_surf) result(ierr)
#endif
        use num_vars, only: no_plots
        use grid_vars, only: dealloc_grid
        use grid_ops, only: trim_grid
#if ldebug
        use num_vars, only: norm_disc_prec_sol, use_pol_flux_F
        use utilities, only: calc_deriv, con2dis, c
#endif
        
        character(*), parameter :: rout_name = 'plot_X_vec'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid
        type(grid_type), intent(in) :: grid_X                                   ! perturbation grid
        type(grid_type), intent(in) :: grid_sol                                 ! solution grid
        type(eq_type), intent(in) :: eq                                         ! equilibrium variables
#if ldebug
        type(met_type), intent(in) :: met                                       ! metric variables
#endif
        type(X_1_type), intent(in) :: X                                         ! perturbation variables
        type(sol_type), intent(in) :: sol                                       ! solution variables
        real(dp), intent(in) :: XYZ(:,:,:,:)                                    ! X, Y and Z of extended eq_grid
        integer, intent(in) :: X_id                                             ! nr. of Eigenvalue (for output name)
        real(dp), intent(in) :: res_surf(:,:)                                   ! resonant surfaces
        
        ! local variables
        type(grid_type) :: grid_sol_trim                                        ! trimmed sol grid
        integer :: norm_id(2)                                                   ! untrimmed normal indices for trimmed grids
        integer :: kd                                                           ! counter
        integer :: n_t(2)                                                       ! nr. of time steps in quarter period, nr. of quarter periods
        integer :: plot_dim(4)                                                  ! dimensions of plot
        integer :: plot_offset(4)                                               ! local offset of plot
        integer :: col                                                          ! collection type for HDF5 plot
        real(dp), allocatable :: time(:)                                        ! fraction of Alfvén time
        real(dp), allocatable :: XYZ_plot(:,:,:,:,:)                            ! copies of XYZ
        complex(dp) :: omega                                                    ! sqrt of Eigenvalue
        complex(dp), allocatable :: f_plot(:,:,:,:,:)                           ! the function to plot
        character(len=max_str_ln) :: err_msg                                    ! error message
        character(len=max_str_ln) :: var_name(2)                                ! name of variable that is plot
        character(len=max_str_ln) :: file_name(2)                               ! name of file
        character(len=max_str_ln) :: description(2)                             ! description
#if ldebug
        real(dp), allocatable :: loc_r_eq(:)                                    ! loc_r_F of sol grid interpolated in eq grid
        integer :: nm                                                           ! n (pol. flux) or m (tor. flux)
        integer :: i_lo, i_hi                                                   ! upper and lower index for interpolation of eq grid to sol grid
        integer :: id, jd, ld                                                   ! counters
        complex(dp), allocatable :: U_inf(:,:,:,:)                              ! ideal ballooning U
        complex(dp), allocatable :: U_inf_prop(:,:,:)                           ! proportional part of U_inf
#endif
        
        ! initialize ierr
        ierr = 0
        
        ! bypass plots if no_plots
        if (no_plots) return
        
        ! tests
        if (size(XYZ,4).ne.3) then
            ierr = 1
            err_msg = 'X, Y and Z needed'
            CHCKERR(err_msg)
        end if
        if (grid_eq%n(1).ne.size(XYZ,1) .or. &
            &grid_eq%n(2).ne.size(XYZ,2) .or. &
            &grid_sol%loc_n_r.ne.size(XYZ,3)) then
            ierr = 1
            err_msg = 'XYZ needs to have the correct dimensions'
            CHCKERR(err_msg)
        end if
        
        ! trim sol grid
        ierr = trim_grid(grid_sol,grid_sol_trim,norm_id)
        CHCKERR('')
        
        ! user output
        call writo('Plotting individual harmonics')
        call lvl_ud(1)
        
        ! plot information about harmonics
        ierr = plot_harmonics(grid_sol_trim,sol,X_id,res_surf)
        CHCKERR('')
        
        call lvl_ud(-1)
        
        ! user output
        call writo('Plotting normal components')
        call lvl_ud(1)
        
        ! set up n_t
        ! if the  Eigenvalue is negative,  the Eigenfunction explodes,  so limit
        ! n_t(2) to 1. If it is positive, the Eigenfunction oscilates, so choose
        ! n_t(2) = 8 for 2 whole periods
        if (realpart(sol%val(X_id)).lt.0) then                                  ! exploding, unstable
            n_t(1) = 1                                                          ! 1 point per quarter period
            n_t(2) = 1                                                          ! 1 quarter period
        else                                                                    ! oscillating, stable
            n_t(1) = 2                                                          ! 2 points per quarter period
            n_t(2) = 4                                                          ! 4 quarter periods
        end if
        
        ! set collection type
        if (product(n_t).gt.1) then
            col = 2                                                             ! temporal collection
        else
            col = 1                                                             ! no collection
        end if
        
        ! set up plot dimensions and local dimensions
        ! Note: The  angular size is taken  from trimmed eq grid but  the normal
        ! size from the trimmed sol grid.
        plot_dim = [grid_eq%n(1), grid_eq%n(2),grid_sol_trim%n(3),product(n_t)]
        plot_offset = [0,0,grid_sol_trim%i_min-1,0]
        !plot_offset = [0,0,grid_sol_trim%i_min-1,product(n_t)]                  !!! THIS WAS WRITTEN FIRST: IS IT CORRECT?!
        
        ! set up copies of XYZ for plot
        allocate(XYZ_plot(grid_eq%n(1),grid_eq%n(2),grid_sol_trim%loc_n_r,&
            &product(n_t),3))
        
        ! For each time step, calculate the time (as fraction of Alfvén
        ! time) and make a copy of XYZ for X, Y and Z
        allocate(time(product(n_t)))
        do kd = 1,product(n_t)
            if (n_t(1).eq.1) then
                time(kd) = (kd-1) * 0.25
            else
                time(kd) = (mod(kd-1,n_t(1))*1._dp/n_t(1) + &
                    &(kd-1)/n_t(1)) * 0.25
            end if
            XYZ_plot(:,:,:,kd,:) = XYZ(:,:,norm_id(1):norm_id(2),:)
        end do
        
        ! calculate omega  = sqrt(X_val)  and make sure  to select  the decaying
        ! solution
        omega = sqrt(sol%val(X_id))
        if (imagpart(omega).gt.0) omega = - omega                               ! exploding solution, not the decaying one
        
        ! calculate  the function  to  plot: normal  and  geodesic component  of
        ! perturbation X_F
        allocate(f_plot(grid_eq%n(1),grid_eq%n(2),grid_sol%loc_n_r,&
            &product(n_t),2))
        ierr = calc_XUQ(grid_eq,grid_X,grid_sol,eq,X,sol,X_id,1,time,&
            &f_plot(:,:,:,:,1))
        CHCKERR('')
        ierr = calc_XUQ(grid_eq,grid_X,grid_sol,eq,X,sol,X_id,2,time,&
            &f_plot(:,:,:,:,2))
        CHCKERR('')
        
#if ldebug
        if (debug_plot_X_vec) then
            call writo('Checking how well U approximates the pure ballooning &
                &result')
            call lvl_ud(1)
            
            ! set up ideal U
            allocate(U_inf(grid_eq%n(1),grid_eq%n(2),grid_sol%loc_n_r,&
                &product(n_t)))
            
            ! get the normal interpolation factors
            allocate(loc_r_eq(grid_sol%loc_n_r))
            do kd = 1,grid_sol%loc_n_r
                ierr = con2dis(grid_sol%loc_r_F(kd),loc_r_eq(kd),&
                    &grid_eq%loc_r_F)
                CHCKERR('')
            end do
            
            ! derive the X vector
            do ld = 1,product(n_t)
                do jd = 1,grid_eq%n(2)
                    do id = 1,grid_eq%n(1)
                        ierr = calc_deriv(f_plot(id,jd,:,ld,1),&
                            &U_inf(id,jd,:,ld),grid_sol%loc_r_F,1,&
                            &norm_disc_prec_sol)
                        CHCKERR('')
                    end do
                end do
            end do
            
            ! set up dummy variable Theta^alpha + q' theta and nm
            allocate(U_inf_prop(grid_eq%n(1),grid_eq%n(2),grid_eq%loc_n_r))
            if (use_pol_flux_F) then
                do kd = 1,grid_eq%loc_n_r
                    U_inf_prop(:,:,kd) = &
                        &met%h_FD(:,:,kd,c([1,2],.true.),0,0,0)/&
                        &met%h_FD(:,:,kd,c([2,2],.true.),0,0,0) + &
                        &eq%q_saf_FD(kd,1)*grid_eq%theta_F(:,:,kd)
                end do
                nm = X%n(1)
            else
                do kd = 1,grid_eq%loc_n_r
                    U_inf_prop(:,:,kd) = &
                        &met%h_FD(:,:,kd,c([1,2],.true.),0,0,0)/&
                        &met%h_FD(:,:,kd,c([2,2],.true.),0,0,0) - &
                        &eq%rot_t_FD(kd,1)*grid_eq%zeta_F(:,:,kd)
                end do
                nm = X%m(1)
            end if
            
            ! multiply by i/n or i/m and add the proportional part
            do ld = 1,product(n_t)
                do kd = 1,grid_sol%loc_n_r
                    i_lo = floor(loc_r_eq(kd))
                    i_hi = ceiling(loc_r_eq(kd))
                    
                    U_inf(:,:,kd,ld) = iu/X%n(1) * U_inf(:,:,kd,ld) - &
                        &f_plot(:,:,kd,ld,1) * (&
                        &U_inf_prop(:,:,i_lo)+(loc_r_eq(kd)-i_lo)*&
                        &(U_inf_prop(:,:,i_hi)-U_inf_prop(:,:,i_lo)))
                end do
            end do
            deallocate(U_inf_prop)
            
            ! plotting real part
            call plot_HDF5('RE X','TEST_RE_X_'//&
                &trim(i2str(X_id)),realpart(f_plot(:,:,:,1,1)),&
                &[grid_eq%n(1:2),grid_sol%n(3)],[0,0,grid_sol%i_min-1])
            call plot_diff_HDF5(realpart(f_plot(:,:,:,1,2)),&
                &realpart(U_inf(:,:,:,1)),'TEST_RE_U_inf_'//trim(i2str(X_id)),&
                &[grid_eq%n(1:2),grid_sol%n(3)],[0,0,grid_sol%i_min-1],&
                &description='To test whether U approximates the ideal &
                &ballooning result',output_message=.true.)
            
            ! plotting imaginary part
            call plot_HDF5('IM X','TEST_IM_X_'//&
                &trim(i2str(X_id)),imagpart(f_plot(:,:,:,1,1)),&
                &[grid_eq%n(1:2),grid_sol%n(3)],[0,0,grid_sol%i_min-1])
            call plot_diff_HDF5(imagpart(f_plot(:,:,:,1,2)),&
                &imagpart(U_inf(:,:,:,1)),'TEST_IM_U_inf_'//trim(i2str(X_id)),&
                &[grid_eq%n(1:2),grid_sol%n(3)],[0,0,grid_sol%i_min-1],&
                &description='To test whether U approximates the ideal &
                &ballooning result',output_message=.true.)
            
            call writo('Checking done')
            call lvl_ud(-1)
        end if
#endif
        
        ! set up var_name, file_name and description
        var_name(1) = 'Normal comp. of EV'
        file_name(1) = trim(i2str(X_id))//'_X_vec'
        description(1) = 'Normal component of solution vector X_vec &
            &for Eigenvalue '//trim(i2str(X_id))//' with omega = '//&
            &trim(r2str(realpart(omega)))
        var_name(2) = 'Geodesic comp. of EV'
        file_name(2) = trim(i2str(X_id))//'_U_vec'
        description(2) = 'Geodesic compoment of solution vector U_vec &
            &for Eigenvalue '//trim(i2str(X_id))//' with omega = '//&
            &trim(r2str(realpart(omega)))
        
        do kd = 1,2
            call plot_HDF5([var_name(kd)],file_name(kd),&
                &realpart(f_plot(:,:,norm_id(1):norm_id(2),:,kd)),&
                &tot_dim=plot_dim,loc_offset=plot_offset,X=XYZ_plot(:,:,:,:,1),&
                &Y=XYZ_plot(:,:,:,:,2),Z=XYZ_plot(:,:,:,:,3),col=col,&
                &description=description(kd))
        end do
        
        ! clean up
        call dealloc_grid(grid_sol_trim)
        
        call lvl_ud(-1)
    contains
        ! plots the harmonics and their maximum in 2D
        ! Note: This routine needs a trimmed grid!
        integer function plot_harmonics(grid_sol,sol,X_id,res_surf) result(ierr)
            use MPI_utilities, only: wait_MPI, get_ser_var
            use output_ops, only: merge_GP
            use num_vars, only: use_pol_flux_F, GP_max_size, rank, no_plots
            use eq_vars, only: max_flux_p_F, max_flux_t_F
            
            character(*), parameter :: rout_name = 'plot_harmonics'
            
            ! input / output
            type(grid_type), intent(in) :: grid_sol                             ! solution grid
            type(sol_type), intent(in) :: sol                                   ! solution variables
            integer, intent(in) :: X_id                                         ! nr. of Eigenvalue (for output name)
            real(dp), intent(in) :: res_surf(:,:)                               ! resonant surfaces
            
            ! local variables
            integer :: id, kd                                                   ! counters
            character(len=max_str_ln) :: file_name                              ! name of file of plots of this proc.
            character(len=max_str_ln) :: plot_title                             ! title for plots
            real(dp) :: norm_factor                                             ! conversion factor max_flux/2pi from flux to normal coordinate
            real(dp), allocatable :: x_plot(:,:)                                ! x values of plot
            real(dp), allocatable :: y_plot(:,:)                                ! y values of plot
            real(dp), allocatable :: X_vec_max(:)                               ! maximum position index of X_vec of rank
            real(dp), allocatable :: res_surf_loc(:,:)                          ! local copy of res_surf
            complex(dp), allocatable :: X_vec_ser(:,:)                          ! serial MPI Eigenvector
            complex(dp), allocatable :: X_vec_ser_loc(:)                        ! local X_vec_ser
            
            ! initialize ierr
            ierr = 0
            
            ! bypass plots if no_plots
            if (no_plots) return
            
            ! set up serial X_vec on master
            if (rank.eq.0) allocate(X_vec_ser(1:sol%n_mod,1:grid_sol%n(3)))
            do id = 1,sol%n_mod
                ierr = get_ser_var(sol%vec(id,norm_id(1):norm_id(2),X_id),&
                    &X_vec_ser_loc,scatter=.true.)
                CHCKERR('')
                if (rank.eq.0) X_vec_ser(id,:) = X_vec_ser_loc
                deallocate(X_vec_ser_loc)
            end do
            
            ! the rest is done only by master
            if (rank.eq.0) then
                ! set up x_plot
                allocate(x_plot(grid_sol%n(3),sol%n_mod))
                do kd = 1,sol%n_mod
                    x_plot(:,kd) = grid_sol%r_F
                end do
                
                ! set up maximum of each mode at midplane
                allocate(X_vec_max(sol%n_mod))
                X_vec_max = 0.0_dp
                do kd = 1,sol%n_mod
                    X_vec_max(kd) = &
                        &grid_sol%r_F(maxloc(abs(realpart(X_vec_ser(kd,:))),1))
                end do
                
                ! scale x_plot and X_vec_max by max_flux/2pi
                if (use_pol_flux_F) then
                    norm_factor = max_flux_p_F/(2*pi)
                else
                    norm_factor = max_flux_t_F/(2*pi)
                end if
                x_plot = x_plot/norm_factor
                X_vec_max = X_vec_max/norm_factor
                
                ! set up file name of this rank and plot title
                file_name = trim(i2str(X_id))//'_EV_midplane_RE'
                plot_title = 'EV - midplane'
                
                ! print real amplitude of harmonics of eigenvector at midplane
                call print_GP_2D(plot_title,file_name,&
                    &realpart(transpose(X_vec_ser)),x=x_plot,draw=.false.)
                
                ! plot in file
                call draw_GP(plot_title,file_name,file_name,sol%n_mod,1,.false.)
                
                ! plot in file using decoupled 3D in GNUPlot if not too big
                if (sol%n_mod*grid_sol%n(3).le.GP_max_size) then
                    call draw_GP(trim(plot_title)//' - 3D',file_name,&
                        &trim(file_name)//'_3D',sol%n_mod,3,.false.)
                end if
                
                ! plot using HDF5
                call plot_HDF5(trim(plot_title),trim(file_name),&
                    &reshape(realpart(transpose(X_vec_ser)),&
                    &[1,grid_sol%n(3),sol%n_mod]),y=reshape(x_plot,&
                    &[1,grid_sol%n(3),sol%n_mod]))
                
                ! set up file name of this rank and plot title
                file_name = trim(i2str(X_id))//'_EV_midplane_IM'
                plot_title = 'EV - midplane'
                
                ! print imag amplitude of harmonics of eigenvector at midplane
                call print_GP_2D(plot_title,file_name,&
                    &imagpart(transpose(X_vec_ser)),x=x_plot,draw=.false.)
                
                ! plot in file
                call draw_GP(plot_title,file_name,file_name,sol%n_mod,1,.false.)
                
                ! plot in file using decoupled 3D in GNUPlot if not too big
                if (sol%n_mod*grid_sol%n(3).le.GP_max_size) then
                    call draw_GP(trim(plot_title)//' - 3D',file_name,&
                        &trim(file_name)//'_3D',sol%n_mod,3,.false.)
                end if
                
                ! plot using HDF5
                call plot_HDF5(trim(plot_title),trim(file_name),&
                    &reshape(imagpart(transpose(X_vec_ser)),&
                    &[1,grid_sol%n(3),sol%n_mod]),y=reshape(x_plot,&
                    &[1,grid_sol%n(3),sol%n_mod]))
                
                ! deallocate variables
                deallocate(x_plot)
                
                ! only plot maximum if resonant surfaces found
                if (size(res_surf,1).gt.0) then
                    ! set up file name of this rank and plot title
                    file_name = trim(i2str(X_id))//'_EV_max'
                    plot_title = trim(i2str(X_id))//' EV - maximum of modes'
                    
                    ! set up local res_surf
                    allocate(res_surf_loc(size(res_surf,1),size(res_surf,2)))
                    res_surf_loc = res_surf
                    
                    ! scale res_surf_loc(2) by max_flux/2pi
                    res_surf_loc(:,2) = res_surf_loc(:,2)/norm_factor
                    
                    ! set up x_plot and y_plot
                    allocate(x_plot(sol%n_mod,2))
                    allocate(y_plot(sol%n_mod,2))
                    x_plot(:,1) = X_vec_max
                    y_plot(:,1) = [(kd*1._dp,kd=1,sol%n_mod)]
                    x_plot(1:size(res_surf_loc,1),2) = res_surf_loc(:,2)
                    x_plot(size(res_surf_loc,1)+1:sol%n_mod,2) = &
                        &res_surf_loc(size(res_surf_loc,1),2)
                    y_plot(1:size(res_surf_loc,1),2) = res_surf_loc(:,1)
                    y_plot(size(res_surf_loc,1)+1:sol%n_mod,2) = &
                        &res_surf_loc(size(res_surf_loc,1),1)
                    
                    ! plot the maximum at midplane
                    call print_GP_2D(plot_title,file_name,y_plot,x=x_plot,&
                        &draw=.false.)
                    
                    ! draw plot in file
                    call draw_GP(plot_title,file_name,file_name,2,1,.false.,&
                        &extra_ops='set xrange ['//&
                        &trim(r2str(grid_sol%r_F(1)/norm_factor))//':'//&
                        &trim(r2str(grid_sol%r_F(grid_sol%n(3))/norm_factor))//&
                        &']')
                end if
            end if
            
            ! synchronize processes
            ierr = wait_MPI()
            CHCKERR('')
        end function plot_harmonics
    end function plot_X_vec
    
    ! Decomposes  the  plasma potential  and  kinetic energy  in its  individual
    ! terms for an individual Eigenvalue.
    ! Use  is made of  variables representing  the potential and  kinetic energy
    ! [ADD REF]:
    !   - E_pot:
    !       + normal line bending term: 1/mu_0 Q_n^2/h22,
    !       + geodesic line bending term: 1/mu_0 J^2 h22/g33 Q_g^2,
    !       + normal ballooning term: -2 p' X^2 kappa_n,
    !       + geodesic ballooning term: -2 p' X U* kappa_g,
    !       + normal kink term: -sigma X*Q_g,
    !       + geodesic kink term: sigma U*Q_n,
    !   - E_kin:
    !       + normal kinetic term: rho X^2/h22,
    !       + geodesic kinetic term: rho J^2 h22/g33 U^2.
    ! The energy  terms are calculated  normally on the sol  grid, interpolating
    ! the quantities defined on the eq grid, and angularly in the eq grid.
    ! Optionally,  the  results can  be  plotted by  providing  X, Y  and Z.  By
    ! default, they are instead written to an output file.
    ! Also, the  fraction between potential and kinetic energy  can be returned,
    ! compared with the Eigenvalue.
    integer function decompose_energy(grid_eq,grid_X,grid_sol,eq,met,X,sol,&
        &X_id,log_i,B_aligned,XYZ,X_val_comp) result(ierr)
        use grid_vars, only: dealloc_grid
        use grid_ops, only: trim_grid
        use num_vars, only: rank, no_plots
        
        character(*), parameter :: rout_name = 'decompose_energy'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid
        type(grid_type), intent(in) :: grid_X                                   ! perturbation grid
        type(grid_type), intent(in) :: grid_sol                                 ! solution grid
        type(eq_type), intent(in) :: eq                                         ! equilibrium variables
        type(met_type), intent(in) :: met                                       ! metric variables
        type(X_1_type), intent(in) :: X                                         ! perturbation variables
        type(sol_type), intent(in) :: sol                                       ! solution variables
        integer, intent(in) :: X_id                                             ! nr. of Eigenvalue
        integer, intent(in) :: log_i                                            ! file number of log file
        logical, intent(in) :: B_aligned                                        ! whether grid is field-aligned
        real(dp), intent(in), optional :: XYZ(:,:,:,:)                          ! X, Y and Z for plotting
        complex(dp), intent(inout), optional :: X_val_comp(2,2)                 ! comparison of EV and energy fraction
        
        ! local variables
        integer :: norm_id(2)                                                   ! untrimmed normal indices for trimmed grids
        integer :: kd                                                           ! counter
        integer :: loc_dim(3)                                                   ! local dimension
        integer :: tot_dim(3)                                                   ! total dimensions
        integer :: loc_offset(3)                                                ! local offsets
        type(grid_type) :: grid_sol_trim                                        ! trimmed sol grid
        complex(dp), allocatable, target :: E_pot(:,:,:,:)                      ! potential energy
        complex(dp), allocatable, target :: E_kin(:,:,:,:)                      ! kinetic energy
        complex(dp), allocatable :: E_pot_int(:)                                ! integrated potential energy
        complex(dp), allocatable :: E_kin_int(:)                                ! integrated kinetic energy
        complex(dp), pointer :: E_pot_trim(:,:,:,:) => null()                   ! trimmed part of E_pot
        complex(dp), pointer :: E_kin_trim(:,:,:,:) => null()                   ! trimmed part of E_kin
        real(dp), allocatable, target :: X_tot(:,:,:,:), Y_tot(:,:,:,:), &
            &Z_tot(:,:,:,:)                                                     ! multiple copies of X, Y and Z for collection plot
        real(dp), pointer :: X_tot_trim(:,:,:,:) => null()                      ! trimmed part of X
        real(dp), pointer :: Y_tot_trim(:,:,:,:) => null()                      ! trimmed part of Y
        real(dp), pointer :: Z_tot_trim(:,:,:,:) => null()                      ! trimmed part of Z
        character(len=max_str_ln), allocatable :: var_names_pot(:)              ! name of potential energy variables
        character(len=max_str_ln), allocatable :: var_names_kin(:)              ! name of kinetic energy variables
        character(len=max_str_ln), allocatable :: var_names(:)                  ! name of other variables that are plot
        character(len=max_str_ln) :: file_name                                  ! name of file
        character(len=max_str_ln) :: description                                ! description
        character(len=max_str_ln) :: format_val                                 ! format
        
        ! initialize ierr
        ierr = 0
        
        ! user output
        call writo('Calculate energy terms')
        call lvl_ud(1)
        
        ierr = calc_E(grid_eq,grid_X,grid_sol,eq,met,X,sol,B_aligned,X_id,&
            &E_pot,E_kin,E_pot_int,E_kin_int)
        CHCKERR('')
        
        call lvl_ud(-1)
        
        ! set X_val_comp if wanted
        if (present(X_val_comp)) then
            X_val_comp(:,1) = complex(X_id*1._dp,0)
            X_val_comp(:,2) = [sol%val(X_id),sum(E_pot_int)/sum(E_kin_int)]
        end if
        
        if (.not.present(XYZ)) then                                             ! writing to log file
            ! user output
            call writo('Write to log file')
            call lvl_ud(1)
            
            ! only master
            if (rank.eq.0) then
                ! set up format string:
                !   row 1: EV, E_pot/E_kin
                !   row 2: E_pot, E_kin
                !   row 3: E_kin(1), E_kin(2)
                !   row 4: E_pot(1), E_pot(2)
                !   row 5: E_pot(3), E_pot(4)
                !   row 6: E_pot(5), E_pot(6)
                format_val = '("  ",ES23.16," ",ES23.16," ",&
                    &ES23.16," ",ES23.16," ",ES23.16," ",ES23.16," ",ES23.16)'
                
                ! write header
                write(log_i,'(A)') '# Eigenvalue '//trim(i2str(X_id))
                
                ! write values
                write(log_i,format_val) &
                    &realpart(sol%val(X_id)),&
                    &realpart(sum(E_pot_int)/sum(E_kin_int)),&
                    &realpart(sum(E_pot_int)),&
                    &realpart(sum(E_kin_int))
                write(log_i,format_val) &
                    &imagpart(sol%val(X_id)),&
                    &imagpart(sum(E_pot_int)/sum(E_kin_int)), &
                    &imagpart(sum(E_pot_int)),&
                    &imagpart(sum(E_kin_int))
                write(log_i,format_val) &
                    &realpart(E_kin_int(1)),&
                    &realpart(E_kin_int(2))
                write(log_i,format_val) &
                    &imagpart(E_kin_int(1)),&
                    &imagpart(E_kin_int(2))
                write(log_i,format_val) &
                    &realpart(E_pot_int(1)),&
                    &realpart(E_pot_int(2)),&
                    &realpart(E_pot_int(3)),&
                    &realpart(E_pot_int(4)),&
                    &realpart(E_pot_int(5)),&
                    &realpart(E_pot_int(6))
                write(log_i,format_val) &
                    &imagpart(E_pot_int(1)),&
                    &imagpart(E_pot_int(2)),&
                    &imagpart(E_pot_int(3)),&
                    &imagpart(E_pot_int(4)),&
                    &imagpart(E_pot_int(5)),&
                    &imagpart(E_pot_int(6))
            end if
            
            call lvl_ud(-1)
        else                                                                    ! plotting
            ! bypass plots if no_plots
            if (no_plots) return
            
            ! user output
            call writo('Preparing plots')
            call lvl_ud(1)
            
            ! set up potential energy variable names
            allocate(var_names_pot(6))
            var_names_pot(1) = '1. normal line bending term ~ Q_n^2'
            var_names_pot(2) = '2. geodesic line bending term ~ Q_g^2'
            var_names_pot(3) = '3. normal ballooning term ~ -p'' X^2 kappa_n'
            var_names_pot(4) = '4. geodesic ballooning term ~ -p'' X U* kappa_g'
            var_names_pot(5) = '5. normal kink term ~ -sigma X*Q_g'
            var_names_pot(6) = '6. geodesic kink term ~ sigma U*Q_n'
            
            ! set up kinetic energy variable names
            allocate(var_names_kin(2))
            var_names_kin(1) = '1. normal kinetic term ~ rho X^2'
            var_names_kin(2) = '2. geodesic kinetic term ~ rho U^2'
            
            ! trim sol grid
            ierr = trim_grid(grid_sol,grid_sol_trim,norm_id)
            CHCKERR('')
            
            ! set plot variables
            loc_dim = [grid_eq%n(1),grid_eq%n(2),grid_sol_trim%loc_n_r]
            tot_dim = [grid_eq%n(1),grid_eq%n(2),grid_sol_trim%n(3)]
            loc_offset = [0,0,grid_sol_trim%i_min-1]
            allocate(X_tot(loc_dim(1),loc_dim(2),loc_dim(3),6))
            allocate(Y_tot(loc_dim(1),loc_dim(2),loc_dim(3),6))
            allocate(Z_tot(loc_dim(1),loc_dim(2),loc_dim(3),6))
            do kd = 1,6
                X_tot(:,:,:,kd) = XYZ(:,:,norm_id(1):norm_id(2),1)
                Y_tot(:,:,:,kd) = XYZ(:,:,norm_id(1):norm_id(2),2)
                Z_tot(:,:,:,kd) = XYZ(:,:,norm_id(1):norm_id(2),3)
            end do
            allocate(var_names(2))
            
            ! point to the trimmed versions of energy
            E_pot_trim => E_pot(:,:,norm_id(1):norm_id(2),:)
            E_kin_trim => E_kin(:,:,norm_id(1):norm_id(2),:)
            
            call lvl_ud(-1)
            
            ! user output
            call writo('Plot energy terms')
            call lvl_ud(1)
            
            ! E_kin
            file_name = trim(i2str(X_id))//'_E_kin'
            description = 'Kinetic energy constituents'
            X_tot_trim => X_tot(:,:,:,1:2)
            Y_tot_trim => Y_tot(:,:,:,1:2)
            Z_tot_trim => Z_tot(:,:,:,1:2)
            call plot_HDF5(var_names_kin,trim(file_name)//'_RE',&
                &realpart(E_kin_trim),tot_dim=[tot_dim,2],&
                &loc_offset=[loc_offset,0],X=X_tot_trim,Y=Y_tot_trim,&
                &Z=Z_tot_trim,description=description)
            call plot_HDF5(var_names_kin,trim(file_name)//'_IM',&
                &imagpart(E_kin_trim),tot_dim=[tot_dim,2],&
                &loc_offset=[loc_offset,0],X=X_tot_trim,Y=Y_tot_trim,&
                &Z=Z_tot_trim,description=description)
            nullify(X_tot_trim,Y_tot_trim,Z_tot_trim)
            
            ! E_pot
            file_name = trim(i2str(X_id))//'_E_pot'
            description = 'Potential energy constituents'
            X_tot_trim => X_tot(:,:,:,1:6)
            Y_tot_trim => Y_tot(:,:,:,1:6)
            Z_tot_trim => Z_tot(:,:,:,1:6)
            call plot_HDF5(var_names_pot,trim(file_name)//'_RE',&
                &realpart(E_pot_trim),tot_dim=[tot_dim,6],&
                &loc_offset=[loc_offset,0],X=X_tot_trim,Y=Y_tot_trim,&
                &Z=Z_tot_trim,description=description)
            call plot_HDF5(var_names_pot,trim(file_name)//'_IM',&
                &imagpart(E_pot_trim),tot_dim=[tot_dim,6],&
                &loc_offset=[loc_offset,0],X=X_tot_trim,Y=Y_tot_trim,&
                &Z=Z_tot_trim,description=description)
            nullify(X_tot_trim,Y_tot_trim,Z_tot_trim)
            
            ! E_stab
            var_names(1) = '1. stabilizing term'
            var_names(2) = '2. destabilizing term'
            file_name = trim(i2str(X_id))//'_E_stab'
            description = '(de)stabilizing energy'
            X_tot_trim => X_tot(:,:,:,1:2)
            Y_tot_trim => Y_tot(:,:,:,1:2)
            Z_tot_trim => Z_tot(:,:,:,1:2)
            call plot_HDF5(var_names,trim(file_name)//'_RE',realpart(reshape(&
                &[sum(E_pot_trim(:,:,:,1:2),4),sum(E_pot_trim(:,:,:,3:6),4)],&
                &[loc_dim(1),loc_dim(2),loc_dim(3),2])),&
                &tot_dim=[tot_dim,2],loc_offset=[loc_offset,0],&
                &X=X_tot_trim,Y=Y_tot_trim,Z=Z_tot_trim,description=description)
            call plot_HDF5(var_names,trim(file_name)//'_IM',imagpart(reshape(&
                &[sum(E_pot_trim(:,:,:,1:2),4),sum(E_pot_trim(:,:,:,3:6),4)],&
                &[loc_dim(1),loc_dim(2),loc_dim(3),2])),&
                &tot_dim=[tot_dim,2],loc_offset=[loc_offset,0],&
                &X=X_tot_trim,Y=Y_tot_trim,Z=Z_tot_trim,description=description)
            
            ! E
            var_names(1) = '1. potential energy'
            var_names(2) = '2. kinetic energy'
            file_name = trim(i2str(X_id))//'_E'
            description = 'total potential and kinetic energy'
            call plot_HDF5(var_names,trim(file_name)//'_RE',realpart(reshape(&
                &[sum(E_pot_trim,4),sum(E_kin_trim,4)],&
                &[loc_dim(1),loc_dim(2),loc_dim(3),2])),&
                &tot_dim=[tot_dim,2],loc_offset=[loc_offset,0],X=X_tot_trim,&
                &Y=Y_tot_trim,Z=Z_tot_trim,description=description)
            call plot_HDF5(var_names,trim(file_name)//'_IM',imagpart(reshape(&
                &[sum(E_pot_trim,4),sum(E_kin_trim,4)],&
                &[loc_dim(1),loc_dim(2),loc_dim(3),2])),&
                &tot_dim=[tot_dim,2],loc_offset=[loc_offset,0],X=X_tot_trim,&
                &Y=Y_tot_trim,Z=Z_tot_trim,description=description)
            
            call lvl_ud(-1)
            
            ! clean up
            nullify(E_kin_trim,E_pot_trim)
            nullify(X_tot_trim,Y_tot_trim,Z_tot_trim)
            call dealloc_grid(grid_sol_trim)
        end if
    end function decompose_energy
    
    ! Calculate the energy terms in the energy decomposition.
    ! Note: see explanation of routine in "decompose_energy"
    integer function calc_E(grid_eq,grid_X,grid_sol,eq,met,X,sol,B_aligned,&
        &X_id,E_pot,E_kin,E_pot_int,E_kin_int) result(ierr)
        use num_vars, only: use_pol_flux_F, n_procs, norm_style, &
            &norm_disc_prec_sol, rank
        use eq_vars, only: vac_perm
        use utilities, only: c, con2dis
        use grid_ops, only: calc_int_vol, trim_grid, untrim_grid
        use grid_vars, only: alpha, dealloc_grid
        use MPI_utilities, only: get_ser_var, wait_MPI
#if ldebug
        use utilities, only: calc_deriv
#endif
        
        character(*), parameter :: rout_name = 'calc_E'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid
        type(grid_type), intent(in) :: grid_X                                   ! perturbation grid
        type(grid_type), intent(in) :: grid_sol                                 ! and solution grid
        type(eq_type), intent(in) :: eq                                         ! equilibrium variables
        type(met_type), intent(in) :: met                                       ! metric variables
        type(X_1_type), intent(in) :: X                                         ! perturbation variables
        type(sol_type), intent(in) :: sol                                       ! solution variables
        logical, intent(in) :: B_aligned                                        ! whether grid is field-aligned
        integer, intent(in) :: X_id                                             ! nr. of Eigenvalue
        complex(dp), intent(inout), allocatable :: E_pot(:,:,:,:)               ! potential energy
        complex(dp), intent(inout), allocatable :: E_kin(:,:,:,:)               ! kinetic energy
        complex(dp), intent(inout), allocatable :: E_pot_int(:)                 ! integrated potential energy
        complex(dp), intent(inout), allocatable :: E_kin_int(:)                 ! integrated kinetic energy
        
        ! local variables
        integer :: norm_id(2)                                                   ! untrimmed normal indices for trimmed grids
        integer :: jd, kd                                                       ! counter
        integer :: i_lo, i_hi                                                   ! upper and lower index for interpolation of eq grid to sol grid
        integer :: loc_dim(3)                                                   ! local dimension
        type(grid_type) :: grid_sol_trim                                        ! trimmed sol grid
        type(grid_type) :: grid_sol_ghost                                       ! ghosted sol grid
        real(dp) :: loc_r_eq                                                    ! loc_r_F of sol grid interpolated in eq grid
        real(dp), allocatable :: h22(:,:,:), g33(:,:,:), J(:,:,:)               ! interpolated h_FD(2,2), g_FD(3,3) and J_FD
        real(dp), allocatable :: kappa_n(:,:,:), kappa_g(:,:,:)                 ! interpolated kappa_n and kappa_g
        real(dp), allocatable :: sigma(:,:,:)                                   ! interpolated sigma
        real(dp), allocatable :: D2p(:,:,:), rho(:,:,:)                         ! interpolated Dpres_FD, rho
        real(dp), allocatable :: ang_1(:,:,:), ang_2(:,:,:), norm(:)            ! coordinates of grid in which to calculate total energy
        complex(dp), allocatable :: XUQ(:,:,:,:)                                ! X, U, Q_n and Q_g
        complex(dp), allocatable :: E_int_tot(:)                                ! integrated potential or kinetic energy of all processes
        character(len=max_str_ln) :: err_msg                                    ! error message
#if ldebug
        real(dp), allocatable :: S(:,:,:)                                       ! interpolated S
        complex(dp), allocatable :: DU(:,:,:)                                   ! D_par U
        real(dp), allocatable :: DU_ALT(:,:,:)                                  ! DU calculated from U
#endif
        
        ! set loc_dim
        loc_dim = [grid_X%n(1:2),grid_X%loc_n_r]                                ! includes ghost regions of width norm_disc_prec_sol
        
        ! allocate local variables
        allocate(h22(loc_dim(1),loc_dim(2),loc_dim(3)))
        allocate(g33(loc_dim(1),loc_dim(2),loc_dim(3)))
        allocate(J(loc_dim(1),loc_dim(2),loc_dim(3)))
        allocate(kappa_n(loc_dim(1),loc_dim(2),loc_dim(3)))
        allocate(kappa_g(loc_dim(1),loc_dim(2),loc_dim(3)))
        allocate(sigma(loc_dim(1),loc_dim(2),loc_dim(3)))
        allocate(D2p(loc_dim(1),loc_dim(2),loc_dim(3)))
        allocate(rho(loc_dim(1),loc_dim(2),loc_dim(3)))
        allocate(ang_1(loc_dim(1),loc_dim(2),loc_dim(3)))
        allocate(ang_2(loc_dim(1),loc_dim(2),loc_dim(3)))
        allocate(norm(loc_dim(3)))
        allocate(XUQ(loc_dim(1),loc_dim(2),loc_dim(3),4))
        allocate(E_pot(loc_dim(1),loc_dim(2),loc_dim(3),6))
        allocate(E_kin(loc_dim(1),loc_dim(2),loc_dim(3),2))
        allocate(E_pot_int(6))
        allocate(E_kin_int(2))
#if ldebug
        if (debug_calc_E .or. debug_DU) then
            allocate(DU(loc_dim(1),loc_dim(2),loc_dim(3)))
            allocate(S(loc_dim(1),loc_dim(2),loc_dim(3)))
            
            ! calculate D_par U
            ierr = calc_XUQ(grid_eq,grid_X,grid_sol,eq,X,sol,X_id,2,0._dp,DU,&
                &met=met,deriv=.true.)
            CHCKERR('')
        end if
#endif
        
        ! iterate over all normal points in sol grid and interpolate
        do kd = 1,loc_dim(3)
            ierr = con2dis(grid_X%loc_r_F(kd),loc_r_eq,grid_eq%loc_r_F)
            CHCKERR('')
            i_lo = floor(loc_r_eq)
            i_hi = ceiling(loc_r_eq)
            
            h22(:,:,kd) = met%h_FD(:,:,i_lo,c([2,2],.true.),0,0,0)+&
                &(loc_r_eq-i_lo)*&
                &(met%h_FD(:,:,i_hi,c([2,2],.true.),0,0,0)&
                &-met%h_FD(:,:,i_lo,c([2,2],.true.),0,0,0))
            g33(:,:,kd) = met%g_FD(:,:,i_lo,c([3,3],.true.),0,0,0)+&
                &(loc_r_eq-i_lo)*&
                &(met%g_FD(:,:,i_hi,c([3,3],.true.),0,0,0)&
                &-met%g_FD(:,:,i_lo,c([3,3],.true.),0,0,0))
            J(:,:,kd) = met%jac_FD(:,:,i_lo,0,0,0)+(loc_r_eq-i_lo)*&
                &(met%jac_FD(:,:,i_hi,0,0,0)&
                &-met%jac_FD(:,:,i_lo,0,0,0))
            kappa_n(:,:,kd) = eq%kappa_n(:,:,i_lo)+(loc_r_eq-i_lo)*&
                &(eq%kappa_n(:,:,i_hi)-eq%kappa_n(:,:,i_lo))
            kappa_g(:,:,kd) = eq%kappa_g(:,:,i_lo)+(loc_r_eq-i_lo)*&
                &(eq%kappa_g(:,:,i_hi)-eq%kappa_g(:,:,i_lo))
            sigma(:,:,kd) = eq%sigma(:,:,i_lo)+(loc_r_eq-i_lo)*&
                &(eq%sigma(:,:,i_hi)-eq%sigma(:,:,i_lo))
            D2p(:,:,kd) = eq%pres_FD(i_lo,1)+(loc_r_eq-i_lo)*&
                &(eq%pres_FD(i_hi,1)-eq%pres_FD(i_lo,1))
            rho(:,:,kd) = eq%rho(i_lo)+(loc_r_eq-i_lo)*&
                &(eq%rho(i_hi)-eq%rho(i_lo))
#if ldebug
            if (debug_calc_E) then
                S(:,:,kd) = eq%S(:,:,i_lo)+(loc_r_eq-i_lo)*&
                    &(eq%S(:,:,i_hi)-eq%S(:,:,i_lo))
            end if
#endif
        end do
        
        ! set angles and norm
        if (B_aligned) then
            if (use_pol_flux_F) then
                ang_1 = grid_X%theta_F
            else
                ang_1 = grid_X%zeta_F
            end if
            ang_2 = alpha
        else
            ang_1 = grid_X%theta_F
            ang_2 = grid_X%zeta_F
        end if
        norm = grid_X%loc_r_F
        
        ! calculate X, U, Q_n and Q_g
        do kd = 1,4
            ierr = calc_XUQ(grid_eq,grid_X,grid_sol,eq,X,sol,X_id,kd,0._dp,&
                &XUQ(:,:,:,kd),met=met)
            CHCKERR('')
        end do
        
        ! calc kinetic energy
        E_kin(:,:,:,1) = rho/h22*XUQ(:,:,:,1)*conjg(XUQ(:,:,:,1))
        select case (norm_style)
            case (1)                                                            ! normalization of full perpendicular component
                E_kin(:,:,:,2) = &
                    &rho*h22*J**2/g33*XUQ(:,:,:,2)*conjg(XUQ(:,:,:,2))
            case (2)                                                            ! normalization of only normal component
                E_kin(:,:,:,2) = 0._dp
            case default
                err_msg = 'No normalization style associated with '//&
                    &trim(i2str(norm_style))
                ierr = 1
                CHCKERR(err_msg)
        end select
        
        ! calc potential energy
        E_pot(:,:,:,1) = 1._dp/vac_perm*1._dp/h22*XUQ(:,:,:,3)*&
            &conjg(XUQ(:,:,:,3))
        E_pot(:,:,:,2) = 1._dp/vac_perm*h22*J**2/g33*XUQ(:,:,:,4)*&
            &conjg(XUQ(:,:,:,4))
        E_pot(:,:,:,3) = -2*D2p*kappa_n*XUQ(:,:,:,1)*conjg(XUQ(:,:,:,1))
        E_pot(:,:,:,4) = -2*D2p*kappa_g*XUQ(:,:,:,1)*conjg(XUQ(:,:,:,2))
        E_pot(:,:,:,5) = -sigma*XUQ(:,:,:,4)*conjg(XUQ(:,:,:,1))
        E_pot(:,:,:,6) = sigma*XUQ(:,:,:,3)*conjg(XUQ(:,:,:,2))
        
#if ldebug
        if (debug_calc_E) then
            call writo('Testing whether the total potential energy is &
                &equal to the alternative form given in [ADD REF]')
            call lvl_ud(1)
            
            ! alternative formulation for E_pot, always real
            E_pot(:,:,:,3) = (-2*D2p*kappa_n+sigma*S)*&
                &XUQ(:,:,:,1)*conjg(XUQ(:,:,:,1))
            E_pot(:,:,:,4) = -sigma/J*2*realpart(XUQ(:,:,:,1)*conjg(DU))
            E_pot(:,:,:,5:6) = 0._dp
            
            call lvl_ud(-1)
        end if
        
        if (debug_X_norm) then
            call writo('Plotting the squared amplitude of the normal &
                &component')
            call lvl_ud(1)
            
            call plot_HDF5('X_norm', 'TEST_X_norm_POST_'//trim(i2str(X_id)),&
                &realpart(XUQ(:,:,:,1)*conjg(XUQ(:,:,:,1))),&
                &tot_dim=grid_X%n,loc_offset=[0,0,grid_X%i_min-1])
            
            call lvl_ud(-1)
        end if
        
        if (debug_DU .and. B_aligned) then
            call writo('Testing whether DU is indeed the parallel &
                &derivative of U')
            call lvl_ud(1)
            
            allocate(DU_ALT(loc_dim(1),loc_dim(2),loc_dim(3)))
            
            ! real part
            do kd = 1,loc_dim(3)
                do jd = 1,loc_dim(2)
                    ierr = calc_deriv(realpart(XUQ(:,jd,kd,2)),&
                        &DU_ALT(:,jd,kd),ang_1(:,jd,kd),1,&
                        &norm_disc_prec_sol)
                    CHCKERR('')
                end do
            end do
            call plot_HDF5('RE U','TEST_RE_U_'//&
                &trim(i2str(X_id)),realpart(XUQ(:,:,:,2)),grid_X%n,&
                &[0,0,grid_sol%i_min-1])
            call plot_diff_HDF5(realpart(DU),DU_ALT,&
                &'TEST_RE_DU_'//trim(i2str(X_id)),grid_X%n,&
                &[0,0,grid_sol%i_min-1],description='To test whether DU is &
                &parallel derivative of U',output_message=.true.)
            
            ! imaginary part
            do kd = 1,loc_dim(3)
                do jd = 1,loc_dim(2)
                    ierr = calc_deriv(imagpart(XUQ(:,jd,kd,2)),&
                        &DU_ALT(:,jd,kd),ang_1(:,jd,kd),1,&
                        &norm_disc_prec_sol)
                    CHCKERR('')
                end do
            end do
            call plot_HDF5('IM U','TEST_IM_U_'//&
                &trim(i2str(X_id)),imagpart(XUQ(:,:,:,2)),grid_X%n,&
                &[0,0,grid_sol%i_min-1])
            call plot_diff_HDF5(imagpart(DU),DU_ALT,'TEST_IM_DU_'//&
                &trim(i2str(X_id)),grid_X%n,[0,0,grid_sol%i_min-1],&
                &description='To test whether DU is parallel derivative &
                &of U',output_message=.true.)
            
            deallocate(DU_ALT)
            
            call lvl_ud(-1)
        end if
#endif
        
        ! trim sol grid
        ierr = trim_grid(grid_sol,grid_sol_trim,norm_id)
        CHCKERR('')
        
        ! add ghost region of width one to the right of the interval
        ierr = untrim_grid(grid_sol_trim,grid_sol_ghost,1)
        CHCKERR('')
        if (rank.lt.n_procs-1) norm_id(2) = norm_id(2)+1                        ! adjust norm_id as well
        
        ! integrate energy using ghosted variables
        ierr = calc_int_vol(ang_1(:,:,norm_id(1):norm_id(2)),&
            &ang_2(:,:,norm_id(1):norm_id(2)),norm(norm_id(1):norm_id(2)),&
            &J(:,:,norm_id(1):norm_id(2)),&
            &E_kin(:,:,norm_id(1):norm_id(2),:),E_kin_int)
        CHCKERR('')
        ierr = calc_int_vol(ang_1(:,:,norm_id(1):norm_id(2)),&
            &ang_2(:,:,norm_id(1):norm_id(2)),norm(norm_id(1):norm_id(2)),&
            &J(:,:,norm_id(1):norm_id(2)),&
            &E_pot(:,:,norm_id(1):norm_id(2),:),E_pot_int)
        CHCKERR('')
        
        ! bundle all processes
        do kd = 1,2
            ierr = get_ser_var([E_kin_int(kd)],E_int_tot,scatter=.true.)
            CHCKERR('')
            E_kin_int(kd) = sum(E_int_tot)
            deallocate(E_int_tot)
        end do
        do kd = 1,6
            ierr = get_ser_var([E_pot_int(kd)],E_int_tot,scatter=.true.)
            CHCKERR('')
            E_pot_int(kd) = sum(E_int_tot)
            deallocate(E_int_tot)
        end do
        
        ! normalize energies
        E_kin = E_kin
        E_pot = E_pot
        E_kin_int = E_kin_int
        E_pot_int = E_pot_int
        
        ! deallocate variables
        call dealloc_grid(grid_sol_trim)
        call dealloc_grid(grid_sol_ghost)
    end function calc_E
    
    ! Print solution quantities to an output file:
    !   - sol:    val, vec
    integer function print_output_sol(grid,sol) result(ierr)
        use num_vars, only: rank, PB3D_name
        use rich, only: rich_info_short
        use HDF5_ops, only: print_HDF5_arrs
        use HDF5_vars, only: var_1D_type
        use grid_ops, only: trim_grid
        use grid_vars, only: dealloc_grid
        
        character(*), parameter :: rout_name = 'print_output_sol'
        
        ! input / output
        type(grid_type), intent(in) :: grid                                     ! solution grid variables
        type(sol_type), intent(in) :: sol                                       ! solution variables
        
        ! local variables
        integer :: norm_id(2)                                                   ! untrimmed normal indices for trimmed grids
        type(var_1D_type), allocatable, target :: sol_1D(:)                     ! 1D equivalent of eq. variables
        type(var_1D_type), pointer :: sol_1D_loc => null()                      ! local element in sol_1D
        type(grid_type) :: grid_trim                                            ! trimmed solution grid
        integer :: id                                                           ! counter
        
        ! initialize ierr
        ierr = 0
        
        ! only write to output file if at least one solution
        if (size(sol%val).gt.0) then
            ! user output
            call writo('Writing solution variables to output file')
            call lvl_ud(1)
            
            ! user output
            call writo('Preparing variables for writing')
            call lvl_ud(1)
            
            ! trim grids
            ierr = trim_grid(grid,grid_trim,norm_id)
            CHCKERR('')
            
            ! Set up the 1D equivalents of the solution variables
            allocate(sol_1D(4))
            
            id = 1
            
            ! RE_X_val
            sol_1D_loc => sol_1D(id); id = id+1
            sol_1D_loc%var_name = 'RE_X_val'
            allocate(sol_1D_loc%tot_i_min(1),sol_1D_loc%tot_i_max(1))
            allocate(sol_1D_loc%loc_i_min(1),sol_1D_loc%loc_i_max(1))
            sol_1D_loc%tot_i_min = [1]
            sol_1D_loc%tot_i_max = [size(sol%val)]
            sol_1D_loc%loc_i_min = [1]
            if (rank.eq.0) then
                sol_1D_loc%loc_i_max = [size(sol%val)]
                allocate(sol_1D_loc%p(size(sol%val)))
                sol_1D_loc%p = realpart(sol%val)
            else
                sol_1D_loc%loc_i_max = [0]
                allocate(sol_1D_loc%p(0))
            end if
            
            ! IM_X_val
            sol_1D_loc => sol_1D(id); id = id+1
            sol_1D_loc%var_name = 'IM_X_val'
            allocate(sol_1D_loc%tot_i_min(1),sol_1D_loc%tot_i_max(1))
            allocate(sol_1D_loc%loc_i_min(1),sol_1D_loc%loc_i_max(1))
            sol_1D_loc%tot_i_min = [1]
            sol_1D_loc%tot_i_max = [size(sol%val)]
            sol_1D_loc%loc_i_min = [1]
            if (rank.eq.0) then
                sol_1D_loc%loc_i_max = [size(sol%val)]
                allocate(sol_1D_loc%p(size(sol%val)))
                sol_1D_loc%p = imagpart(sol%val)
            else
                sol_1D_loc%loc_i_max = [0]
                allocate(sol_1D_loc%p(0))
            end if
            
            ! RE_X_vec
            sol_1D_loc => sol_1D(id); id = id+1
            sol_1D_loc%var_name = 'RE_X_vec'
            allocate(sol_1D_loc%tot_i_min(3),sol_1D_loc%tot_i_max(3))
            allocate(sol_1D_loc%loc_i_min(3),sol_1D_loc%loc_i_max(3))
            sol_1D_loc%loc_i_min = [1,grid_trim%i_min,1]
            sol_1D_loc%loc_i_max = [sol%n_mod,grid_trim%i_max,size(sol%vec,3)]
            sol_1D_loc%tot_i_min = [1,1,1]
            sol_1D_loc%tot_i_max = [sol%n_mod,grid_trim%n(3),size(sol%vec,3)]
            allocate(sol_1D_loc%p(size(sol%vec(:,norm_id(1):norm_id(2),:))))
            sol_1D_loc%p = &
                &reshape(realpart(sol%vec(:,norm_id(1):norm_id(2),:)),&
                &[size(sol%vec(:,norm_id(1):norm_id(2),:))])
            
            ! IM_X_vec
            sol_1D_loc => sol_1D(id); id = id+1
            sol_1D_loc%var_name = 'IM_X_vec'
            allocate(sol_1D_loc%tot_i_min(3),sol_1D_loc%tot_i_max(3))
            allocate(sol_1D_loc%loc_i_min(3),sol_1D_loc%loc_i_max(3))
            sol_1D_loc%loc_i_min = [1,grid_trim%i_min,1]
            sol_1D_loc%loc_i_max = [sol%n_mod,grid_trim%i_max,size(sol%vec,3)]
            sol_1D_loc%tot_i_min = [1,1,1]
            sol_1D_loc%tot_i_max = [sol%n_mod,grid_trim%n(3),size(sol%vec,3)]
            allocate(sol_1D_loc%p(size(sol%vec(:,norm_id(1):norm_id(2),:))))
            sol_1D_loc%p = &
                &reshape(imagpart(sol%vec(:,norm_id(1):norm_id(2),:)),&
                &[size(sol%vec(:,norm_id(1):norm_id(2),:))])
            
            call lvl_ud(-1)
            
            ! user output
            call writo('Writing using HDF5')
            call lvl_ud(1)
            
            ! write
            ierr = print_HDF5_arrs(sol_1D,PB3D_name,'sol'//&
                &trim(rich_info_short()))
            CHCKERR('')
            
            ! clean up
            call dealloc_grid(grid_trim)
            nullify(sol_1D_loc)
            
            ! user output
            call lvl_ud(-1)
            
            ! user output
            call lvl_ud(-1)
            call writo('Solution variables written to output')
        else
            ! user output
            call writo('No solutions to write')
        end if
    end function print_output_sol
end module sol_ops
