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
    use X_vars, only: X_type

    implicit none
    private
    public calc_XUQ, plot_X_vec, decompose_energy
#if ldebug
    public debug_calc_XUQ_arr
#endif
    
    ! interfaces
    interface calc_XUQ
        module procedure calc_XUQ_arr, calc_XUQ_ind
    end interface
    
    ! global variables
#if ldebug
    logical :: debug_calc_XUQ_arr = .false.                                     ! plot debug information for interp_fun_0D_real
#endif

contains
    ! calculates the  normal or geodesic  component, or optionally  the parallel
    ! derivative of the plasma perturbation  or the normal or geodesic component
    ! of the magnetic field perturbation in the perturbation grid for a range in
    ! normalized  time (1  corresponds to  one period).  The variable  XUQ_style
    ! determines which one:
    !   - XUQ_style = 1: X (supports parallel derivative)
    !   - XUQ_style = 2: U (supports parallel derivative)
    !   - XUQ_style = 3: Qn
    !   - XUQ_style = 4: Qg
    ! For Qn  and Qg,  the metric  variables have  to be  provided as  well. The
    ! output is given in the X grid for the normal part. The angular part of the
    ! output is given by  the equilibrium grid. These grids need  to be given in
    ! Flux coordinates.
    ! Note: The angular part of the X grid is neglected as it is assumed to be a
    ! 1D grid.
    integer function calc_XUQ_arr(grid_eq,eq,grid_X,X,X_id,XUQ_style,time,XUQ,&
        &met,deriv) result(ierr)                                                ! (time) array version
        use num_vars, only: use_pol_flux_F, ghost_width_POST
        use utilities, only: con2dis, calc_deriv
#if ldebug
        use utilities, only: calc_int
#endif
        
        character(*), parameter :: rout_name = 'calc_XUQ_arr'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid
        type(eq_type), intent(in) :: eq                                         ! equilibrium variables
        type(grid_type), intent(in) :: grid_X                                   ! perturbation grid
        type(X_type), intent(in) :: X                                           ! perturbation variables
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
        integer :: deriv_ord                                                    ! order of derivative
        complex(dp) :: sqrt_X_val_norm                                          ! normalized sqrt(X_val)
        logical :: deriv_loc                                                    ! local copy of deriv
        complex(dp), allocatable :: DX_vec(:)                                   ! normal derivative of X_vec for a specific mode
        complex(dp), allocatable :: fac_0(:,:,:), fac_1(:,:,:)                  ! factor to be multiplied with X and DX
        complex(dp), allocatable :: par_fac(:)                                  ! multiplicative factor due to parallel derivative
        complex(dp), allocatable :: XUQ_loc(:,:)                                ! complex XUQ without time at a normal point
        real(dp), allocatable :: grp_r_eq(:)                                    ! grp_r_F of X grid interpolated in eq grid
        integer :: i_lo, i_hi                                                   ! upper and lower index for interpolation of eq grid to X grid
#if ldebug
        real(dp), allocatable :: X_vec_ALT(:)                                   ! X_vec calculated from DX_vec
#endif
        
        ! initialize ierr
        ierr = 0
        
        ! set up n_t
        n_t = size(time)
        
        ! tests
        if (grid_eq%n(1).ne.size(XUQ,1) .or. grid_eq%n(2).ne.size(XUQ,2) .or. &
            &grid_X%grp_n_r.ne.size(XUQ,3) .or. n_t.ne.size(XUQ,4)) then
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
        
        ! set up order of derivative
        select case (ghost_width_POST)
            case (:0)                                                           ! no ghost region
                ierr = 1
                err_msg = 'Need ghost region width of at least 1'
                CHCKERR(err_msg)
            case (1)
                deriv_ord = 1                                                   ! 1 point available so derive with order 1
            case (2:)
                deriv_ord = 2                                                   ! 2 point available so derive with order 2
        end select
        
        ! set up local copy of deriv
        deriv_loc = .false.
        if (present(deriv)) deriv_loc = deriv
        
        ! set up normalized sqrt(X_val)
        sqrt_X_val_norm = sqrt(X%val(X_id))
        if (imagpart(sqrt_X_val_norm).gt.0) &
            &sqrt_X_val_norm = - sqrt_X_val_norm                                ! exploding solution, not the decaying one
        sqrt_X_val_norm = sqrt_X_val_norm / abs(sqrt_X_val_norm)                ! normalize
        
        ! calculate  normal  interpolation  tables for  equilibrium  grid  (this
        ! concerns eq variables and met variables)
        allocate(grp_r_eq(grid_X%grp_n_r))
        do kd = 1,grid_X%grp_n_r
            ierr = con2dis(grid_X%grp_r_F(kd),grp_r_eq(kd),grid_eq%grp_r_F)
            CHCKERR('')
        end do
        
        ! Initialize  multiplicative factors fac_0  and fac_1 and  par_fac which
        ! are used in XUQ = fac_0*X + fac_1*DX with fac_0 and fac_1:
        !   1. fac_0 = 1 (nq-m or n-iota m for deriv)
        !      fac_1 = 0
        !   2. fac_0 = U_0 (DU_0 for deriv)
        !      fac_1 = U_1 (DU_1 for deriv)
        !   3. fac_0 = (nq-m)/J or (n-iotam)/J (no deriv)
        !      fac_1 = 0
        !   4. fac_0 = DU_0/J - S (no deriv)
        !      fac_1 = DU_1/J
        ! and par_fac is a helper variable.
        ! Note that if the resolution for X_vec is bad, the numerical derivative
        ! is very  inaccurate, therefore  only smooth (i.e.  physical) solutions
        ! should be looked at.
        allocate(fac_0(grid_eq%n(1),grid_eq%n(2),grid_eq%grp_n_r))
        allocate(fac_1(grid_eq%n(1),grid_eq%n(2),grid_eq%grp_n_r))
        allocate(par_fac(grid_eq%grp_n_r))
        
        ! set up DX_vec
        allocate(DX_vec(grid_X%grp_n_r))
        
        ! initialize XUQ
        XUQ = 0._dp
        
        ! initialize XUQ_loc
        allocate(XUQ_loc(size(XUQ,1),size(XUQ,2)))
        
        ! iterate over all modes
        do jd = 1,X%n_mod
            ! set up parallel multiplicative factor par_fac in normal eq grid
            if (use_pol_flux_F) then
                par_fac = iu*(X%n(jd)*eq%q_saf_FD(:,0)-X%m(jd))
            else
                par_fac = iu*(X%n(jd)-X%m(jd)*eq%rot_t_FD(:,0))
            end if
            
            ! set up multiplicative factors in eq grid
            select case (XUQ_style)
                case (1)                                                        ! calculating X
                    if (deriv_loc) then                                         ! parallel derivative
                        do kd = 1,grid_eq%grp_n_r
                            fac_0(:,:,kd) = par_fac(kd)
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
                        do kd = 1,grid_eq%grp_n_r
                            fac_0(:,:,kd) = par_fac(kd)/met%jac_FD(:,:,kd,0,0,0)
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
                        do kd = 1,grid_eq%grp_n_r
                            fac_0(:,:,kd) = -eq%S(:,:,kd) +&
                                &X%DU_0(:,:,kd,jd)/met%jac_FD(:,:,kd,0,0,0)
                            fac_1(:,:,kd) = X%DU_1(:,:,kd,jd)/&
                                &met%jac_FD(:,:,kd,0,0,0)
                        end do
                    end if
                case default
                    err_msg = 'No XUQ XUQ_style associated with '//&
                        &trim(i2str(XUQ_style))
                    ierr = 1
                    CHCKERR(err_msg)
            end select
            
            ! set up normal derivative of X vec
            if (XUQ_style.eq.2 .or. XUQ_style.eq.4) then
                ierr = calc_deriv(X%vec(jd,:,X_id),DX_vec,&
                    &grid_X%grp_r_F(1:grid_X%grp_n_r),1,deriv_ord)
                CHCKERR('')
#if ldebug
                if (debug_calc_XUQ_arr) then
                    call writo('For mode '//trim(i2str(jd))//', testing &
                        &whether DX_vec is correct by comparing its integral &
                        &with original X_vec')
                    allocate(X_vec_ALT(1:grid_X%grp_n_r))
                    ierr = calc_int(realpart(DX_vec),&
                        &grid_X%grp_r_F(1:grid_X%grp_n_r),X_vec_ALT)
                    CHCKERR('')
                    call print_GP_2D('TEST_RE_DX_vec','',reshape(&
                        &[realpart(X%vec(jd,:,X_id)),X_vec_ALT],&
                        &[grid_X%grp_n_r,2]),x=reshape(&
                        &[grid_X%grp_r_F(1:grid_X%grp_n_r),&
                        &grid_X%grp_r_F(1:grid_X%grp_n_r)],&
                        &[grid_X%grp_n_r,2]))
                    ierr = calc_int(imagpart(DX_vec),&
                        &grid_X%grp_r_F(1:grid_X%grp_n_r),X_vec_ALT)
                    CHCKERR('')
                    call print_GP_2D('TEST_IM_DX_vec','',reshape(&
                        &[imagpart(X%vec(jd,:,X_id)),X_vec_ALT],&
                        &[grid_X%grp_n_r,2]),x=reshape(&
                        &[grid_X%grp_r_F(1:grid_X%grp_n_r),&
                        &grid_X%grp_r_F(1:grid_X%grp_n_r)],&
                        &[grid_X%grp_n_r,2]))
                    deallocate(X_vec_ALT)
                end if
#endif
            else
                DX_vec = 0._dp
            end if
            
            ! iterate over all normal points in X grid (of this group)
            do kd = 1,grid_X%grp_n_r
                ! set lower and upper index for this normal point
                i_lo = floor(grp_r_eq(kd))
                i_hi = ceiling(grp_r_eq(kd))
                
                ! set up loc complex XUQ without time at this normal point
                XUQ_loc(:,:) = exp(iu*(&
                    &X%n(jd)*&
                    &(grid_eq%zeta_F(:,:,i_lo)+(grp_r_eq(kd)-i_lo)*&
                    &(grid_eq%zeta_F(:,:,i_hi)-grid_eq%zeta_F(:,:,i_lo)))&
                    &- X%m(jd)*&
                    &(grid_eq%theta_F(:,:,i_lo)+(grp_r_eq(kd)-i_lo)*&
                    &(grid_eq%theta_F(:,:,i_hi)-grid_eq%theta_F(:,:,i_lo)))&
                    &)) * ( &
                    &X%vec(jd,kd,X_id) * &
                    &(fac_0(:,:,i_lo)+(grp_r_eq(kd)-i_lo)*&
                    &(fac_0(:,:,i_hi)-fac_0(:,:,i_lo))) &
                    &+ DX_vec(kd) * &
                    &(fac_1(:,:,i_lo)+(grp_r_eq(kd)-i_lo)*&
                    &(fac_1(:,:,i_hi)-fac_1(:,:,i_lo))) &
                    &)
                
                ! iterate over time steps
                do id = 1,n_t
                    ! add current mode
                    XUQ(:,:,kd,id) = XUQ(:,:,kd,id) + XUQ_loc * &
                        &exp(iu*sqrt_X_val_norm*time(id)*2*pi)
                end do
            end do
        end do
    end function calc_XUQ_arr
    integer function calc_XUQ_ind(grid_eq,eq,grid_X,X,X_id,XUQ_style,time,XUQ,&
        &met,deriv) result(ierr)                                                ! (time) individual version
        character(*), parameter :: rout_name = 'calc_XUQ_ind'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! equilibirum grid
        type(eq_type), intent(in) :: eq                                         ! equilibrium variables
        type(grid_type), intent(in) :: grid_X                                   ! perturbation grid
        type(X_type), intent(in) :: X                                           ! perturbation variables
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
        ierr = calc_XUQ_arr(grid_eq,eq,grid_X,X,X_id,XUQ_style,[time],&
            &XUQ_arr,met,deriv)
        CHCKERR('')
        
        ! copy array to individual XUQ
        XUQ = XUQ_arr(:,:,:,1)
        
        ! deallocate array version of XUQ
        deallocate(XUQ_arr)
    end function calc_XUQ_ind
    
    ! Plots  Eigenvectors  using  the  angular  part  of  the  the  provided
    ! equilibrium  grid and  the normal  part of  the provided  perturbation
    ! grid.
    integer function plot_X_vec(grid_eq,eq,grid_X,X,XYZ,X_id,res_surf) &
        &result(ierr)
        use output_ops, only: plot_HDF5
        use grid_ops, only: trim_grid
        use grid_vars, only: dealloc_grid
        use num_vars, only: grp_rank, ghost_width_POST
        
        character(*), parameter :: rout_name = 'plot_X_vec'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid
        type(eq_type), intent(in) :: eq                                         ! equilibrium variables
        type(grid_type), intent(in) :: grid_X                                   ! perturbation grid
        type(X_type), intent(in) :: X                                           ! perturbation variables
        real(dp), intent(in) :: XYZ(:,:,:,:)                                    ! X, Y and Z of extended eq_grid
        integer, intent(in) :: X_id                                             ! nr. of Eigenvalue (for output name)
        real(dp), intent(in) :: res_surf(:,:)                                   ! resonant surfaces
        
        ! local variables
        integer :: kd                                                           ! counter
        integer :: n_t(2)                                                       ! nr. of time steps in quarter period, nr. of quarter periods
        type(grid_type) :: grid_X_trim                                          ! trimmed X grid
        character(len=max_str_ln) :: err_msg                                    ! error message
        complex(dp), allocatable :: f_plot(:,:,:,:,:)                           ! the function to plot
        real(dp), allocatable :: time(:)                                        ! fraction of Alfvén time
        character(len=max_str_ln) :: var_name(2)                                ! name of variable that is plot
        character(len=max_str_ln) :: file_name(2)                               ! name of file
        character(len=max_str_ln) :: description(2)                             ! description
        complex(dp) :: omega                                                    ! sqrt of Eigenvalue
        integer :: plot_dim(4)                                                  ! dimensions of plot
        integer :: plot_offset(4)                                               ! group offset of plot
        integer :: col                                                          ! collection type for HDF5 plot
        integer :: norm_ut(2)                                                   ! untrimmed normal indices for trimmed grids
        real(dp), allocatable :: X_plot(:,:,:,:)                                ! copies of XYZ(1)
        real(dp), allocatable :: Y_plot(:,:,:,:)                                ! copies of XYZ(2)
        real(dp), allocatable :: Z_plot(:,:,:,:)                                ! copies of XYZ(3)
        
        ! initialize ierr
        ierr = 0
        
        ! tests
        if (size(XYZ,4).ne.3) then
            ierr = 1
            err_msg = 'X, Y and Z needed'
            CHCKERR(err_msg)
        end if
        if (grid_eq%n(1).ne.size(XYZ,1) .or. &
            &grid_eq%n(2).ne.size(XYZ,2) .or. &
            &grid_X%grp_n_r.ne.size(XYZ,3)) then
            ierr = 1
            err_msg = 'XYZ needs to have the correct dimensions'
            CHCKERR(err_msg)
        end if
        
        ! trim X grid
        ierr = trim_grid(grid_X,grid_X_trim,shift_grid=ghost_width_POST)
        CHCKERR('')
        norm_ut = [1,grid_X_trim%grp_n_r]
        if (grp_rank.gt.0) norm_ut = norm_ut + ghost_width_POST                 ! everything's shifted by ghost_width_POST if not first process
        
        ! user output
        call writo('Plotting individual harmonics')
        call lvl_ud(1)
        
        ! plot information about harmonics
        ierr = plot_harmonics(grid_X_trim,X,X_id,res_surf)
        CHCKERR('')
        
        call lvl_ud(-1)
        
        ! user output
        call writo('Plotting normal components')
        call lvl_ud(1)
        
        ! set up n_t
        ! if the  Eigenvalue is negative,  the Eigenfunction explodes,  so limit
        ! n_t(2) to 1. If it is positive, the Eigenfunction oscilates, so choose
        ! n_t(2) = 8 for 2 whole periods
        if (realpart(X%val(X_id)).lt.0) then                                    ! exploding, unstable
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
        
        ! set up plot dimensions and group dimensions
        ! Note: The  angular size is taken  from trimmed eq grid but  the normal
        ! size from the trimmed X grid.
        plot_dim = [grid_eq%n(1), grid_eq%n(2),grid_X_trim%n(3),product(n_t)]
        plot_offset = [0,0,grid_X_trim%i_min-1,product(n_t)]
        
        ! set up copies  of XYZ(1), XYZ(2) and XYZ(3) in  X, Y and Z
        ! for plot
        allocate(X_plot(grid_eq%n(1),grid_eq%n(2),grid_X_trim%grp_n_r,&
            &product(n_t)))
        allocate(Y_plot(grid_eq%n(1),grid_eq%n(2),grid_X_trim%grp_n_r,&
            &product(n_t)))
        allocate(Z_plot(grid_eq%n(1),grid_eq%n(2),grid_X_trim%grp_n_r,&
            &product(n_t)))
        
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
            X_plot(:,:,:,kd) = XYZ(:,:,norm_ut(1):norm_ut(2),1)
            Y_plot(:,:,:,kd) = XYZ(:,:,norm_ut(1):norm_ut(2),2)
            Z_plot(:,:,:,kd) = XYZ(:,:,norm_ut(1):norm_ut(2),3)
        end do
        
        ! calculate omega  = sqrt(X_val)  and make sure  to select  the decaying
        ! solution
        omega = sqrt(X%val(X_id))
        if (imagpart(omega).gt.0) omega = - omega                               ! exploding solution, not the decaying one
        
        ! calculate  the function  to  plot: normal  and  geodesic component  of
        ! perturbation X_F
        allocate(f_plot(grid_eq%n(1),grid_eq%n(2),grid_X%grp_n_r,&
            &product(n_t),2))
        ierr = calc_XUQ(grid_eq,eq,grid_X,X,X_id,1,time,f_plot(:,:,:,:,1))
        CHCKERR('')
        ierr = calc_XUQ(grid_eq,eq,grid_X,X,X_id,2,time,f_plot(:,:,:,:,2))
        CHCKERR('')
        
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
                &realpart(f_plot(:,:,norm_ut(1):norm_ut(2),:,kd)),&
                &tot_dim=plot_dim,grp_offset=plot_offset,X=X_plot,Y=Y_plot,&
                &Z=Z_plot,col=col,description=description(kd))
        end do
        
        ! deallocate local variables
        deallocate(time)
        deallocate(f_plot)
        deallocate(X_plot,Y_plot,Z_plot)
        call dealloc_grid(grid_X_trim)
        
        call lvl_ud(-1)
    contains
        ! plots the harmonics and their maximum in 2D
        ! Note: This routine needs a trimmed grid!
        integer function plot_harmonics(grid_X,X,X_id,res_surf) result(ierr)
            use MPI_utilities, only: wait_MPI, get_ghost_arr, get_ser_var
            use output_ops, only: merge_GP
            use num_vars, only: plt_rank, use_pol_flux_F, GP_max_size
            use eq_vars, only: max_flux_p_F, max_flux_t_F
            
            character(*), parameter :: rout_name = 'plot_harmonics'
            
            ! input / output
            type(grid_type), intent(in) :: grid_X                               ! perturbation grid
            type(X_type), intent(in) :: X                                       ! perturbation variables
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
            
            ! set up serial X_vec on group master
            if (plt_rank.eq.0) allocate(X_vec_ser(1:X%n_mod,1:grid_X%n(3)))
            do id = 1,X%n_mod
                ierr = get_ser_var(X%vec(id,norm_ut(1):norm_ut(2),X_id),&
                    &X_vec_ser_loc,scatter=.true.)
                CHCKERR('')
                if (plt_rank.eq.0) X_vec_ser(id,:) = X_vec_ser_loc
                deallocate(X_vec_ser_loc)
            end do
            
            ! the rest is done only by plot master
            if (plt_rank.eq.0) then
                ! set up x_plot
                allocate(x_plot(grid_X%n(3),X%n_mod))
                do kd = 1,X%n_mod
                    x_plot(:,kd) = grid_X%r_F
                end do
                
                ! set up local res_surf
                allocate(res_surf_loc(size(res_surf,1),size(res_surf,2)))
                res_surf_loc = res_surf
                
                ! set up maximum of each mode at midplane
                allocate(X_vec_max(X%n_mod))
                X_vec_max = 0.0_dp
                do kd = 1,X%n_mod
                    X_vec_max(kd) = &
                        &grid_X%r_F(maxloc(abs(realpart(X_vec_ser(kd,:))),1))
                end do
                
                ! scale x_plot, X_vec_max and res_surf_loc(2) by max_flux/2pi
                if (use_pol_flux_F) then
                    norm_factor = max_flux_p_F/(2*pi)
                else
                    norm_factor = max_flux_t_F/(2*pi)
                end if
                x_plot = x_plot/norm_factor
                X_vec_max = X_vec_max/norm_factor
                res_surf_loc(:,2) = res_surf_loc(:,2)/norm_factor
                
                ! set up file name of this rank and plot title
                file_name = trim(i2str(X_id))//'_EV_midplane'
                plot_title = 'EV - midplane'
                
                ! print amplitude of harmonics of eigenvector at midplane
                call print_GP_2D(plot_title,file_name,&
                    &realpart(transpose(X_vec_ser)),x=x_plot,draw=.false.)
                
                ! plot in file
                call draw_GP(plot_title,file_name,file_name,X%n_mod,1,.false.)
                
                ! plot in file using decoupled 3D in GNUPlot if not too big
                if (X%n_mod*grid_X%n(3).le.GP_max_size) then
                    call draw_GP(trim(plot_title)//' - 3D',file_name,&
                        &trim(file_name)//'_3D',X%n_mod,3,.false.)
                end if
                
                ! plot using HDF5
                call plot_HDF5(trim(plot_title),trim(file_name),&
                    &reshape(realpart(transpose(X_vec_ser)),&
                    &[1,grid_X%n(3),X%n_mod]),y=reshape(x_plot,&
                    &[1,grid_X%n(3),X%n_mod]))
                
                ! deallocate variables
                deallocate(x_plot)
                
                ! set up file name of this rank and plot title
                file_name = trim(i2str(X_id))//'_EV_max'
                plot_title = trim(i2str(X_id))//' EV - maximum of modes'
                
                ! set up x_plot and y_plot
                allocate(x_plot(X%n_mod,2))
                allocate(y_plot(X%n_mod,2))
                x_plot(:,1) = X_vec_max
                y_plot(:,1) = [(kd*1._dp,kd=1,X%n_mod)]
                x_plot(1:size(res_surf_loc,1),2) = res_surf_loc(:,2)
                x_plot(size(res_surf_loc,1)+1:X%n_mod,2) = &
                    &res_surf_loc(size(res_surf_loc,1),2)
                y_plot(1:size(res_surf_loc,1),2) = res_surf_loc(:,1)
                y_plot(size(res_surf_loc,1)+1:X%n_mod,2) = &
                    &res_surf_loc(size(res_surf_loc,1),1)
                
                ! plot the maximum at midplane
                call print_GP_2D(plot_title,file_name,y_plot,x=x_plot,&
                    &draw=.false.)
                
                ! draw plot in file
                call draw_GP(plot_title,file_name,file_name,2,1,.false.,&
                    &extra_ops='set xrange ['//&
                    &trim(r2str(grid_X%r_F(1)/norm_factor))//':'//&
                    &trim(r2str(grid_X%r_F(grid_X%n(3))/norm_factor))//']')
                
                ! deallocate
                deallocate(x_plot,X_vec_ser)
            end if
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
    ! The energy terms are calculated normally  on the X grid, interpolating the
    ! quantities defined on the eq grid, and angularly in the eq grid.
    ! Optionally, the results can be plotted  by providing PB3D_plot, where X, Y
    ! and Z can be provided as well. If not plotting, they are ignored.
    integer function decompose_energy(PB3D,X_id,log_i,PB3D_plot,XYZ) &
        &result(ierr)
        use PB3D_vars, only: PB3D_type
        use grid_ops, only: trim_grid
        use grid_vars, only: dealloc_grid
        use num_vars, only: ghost_width_POST, grp_rank
        
        character(*), parameter :: rout_name = 'decompose_energy'
        
        ! input / output
        type(PB3D_type), intent(in) :: PB3D                                     ! field-aligned PB3D variables
        integer, intent(in) :: X_id                                             ! nr. of Eigenvalue
        integer, intent(in) :: log_i                                            ! file number of log file
        type(PB3D_type), intent(in), optional :: PB3D_plot                      ! optionally provide plot variables
        real(dp), intent(in), optional :: XYZ(:,:,:,:)                          ! X, Y and Z for plotting
        
        ! local variables
        logical :: plot_en                                                      ! whether plot variables are provided
        integer :: kd                                                           ! counter
        integer :: grp_dim(3)                                                   ! group dimension
        integer :: tot_dim(3)                                                   ! total dimensions
        integer :: grp_offset(3)                                                ! group offsets
        integer :: norm_ut(2)                                                   ! untrimmed normal indices for trimmed grids
        type(grid_type) :: grid_X_trim                                          ! trimmed X grid
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
        
        ! set plot_en
        plot_en = .false.
        if (present(PB3D_plot)) plot_en = .true.
        
        ! user output
        call writo('Prepare calculations')
        call lvl_ud(1)
        
        ! set up format string:
        !   row 1: EV, E_pot/E_kin
        !   row 2: E_pot, E_kin
        !   row 3: E_kin(1), E_kin(2)
        !   row 4: E_pot(1), E_pot(2)
        !   row 5: E_pot(3), E_pot(4)
        !   row 6: E_pot(5), E_pot(6)
        if (grp_rank.eq.0) format_val = '("  ",ES23.16," ",ES23.16," ",&
            &ES23.16," ",ES23.16," ",ES23.16," ",ES23.16," ",ES23.16)'
        
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
        
        call lvl_ud(-1)
        
        call writo('Calculate energy terms')
        call lvl_ud(1)
        
        ierr = calc_E(PB3D,.true.,X_id,E_pot,E_kin,E_pot_int,E_kin_int)
        CHCKERR('')
        
        call lvl_ud(-1)
        
        call writo('Write to log file')
        call lvl_ud(1)
        if (grp_rank.eq.0) write(log_i,'(A)') '# Eigenvalue '//trim(i2str(X_id))
        
        if (grp_rank.eq.0) write(log_i,format_val) &
            &realpart(PB3D%X%val(X_id)),&
            &realpart(sum(E_pot_int)/sum(E_kin_int)),&
            &realpart(sum(E_kin_int)),&
            &realpart(sum(E_pot_int))
        if (grp_rank.eq.0) write(log_i,format_val) &
            &imagpart(PB3D%X%val(X_id)),&
            &imagpart(sum(E_pot_int)/sum(E_kin_int)), &
            &imagpart(sum(E_kin_int)),&
            &imagpart(sum(E_pot_int))
        if (grp_rank.eq.0) write(log_i,format_val) &
            &realpart(E_kin_int(1)),&
            &realpart(E_kin_int(2))
        if (grp_rank.eq.0) write(log_i,format_val) &
            &imagpart(E_kin_int(1)),&
            &imagpart(E_kin_int(2))
        if (grp_rank.eq.0) write(log_i,format_val) &
            &realpart(E_pot_int(1)),&
            &realpart(E_pot_int(2)),&
            &realpart(E_pot_int(3)),&
            &realpart(E_pot_int(4)),&
            &realpart(E_pot_int(5)),&
            &realpart(E_pot_int(6))
        if (grp_rank.eq.0) write(log_i,format_val) &
            &imagpart(E_pot_int(1)),&
            &imagpart(E_pot_int(2)),&
            &imagpart(E_pot_int(3)),&
            &imagpart(E_pot_int(4)),&
            &imagpart(E_pot_int(5)),&
            &imagpart(E_pot_int(6))
        
        call lvl_ud(-1)
        
        ! deallocate variables
        deallocate(E_pot,E_kin,E_pot_int,E_kin_int)
        
        if (plot_en) then
            ! trim X grid
            ierr = trim_grid(PB3D_plot%grid_X,grid_X_trim,&
                &shift_grid=ghost_width_POST)
            CHCKERR('')
            
            ! set up norm_ut
            norm_ut = [1,grid_X_trim%grp_n_r]
            if (grp_rank.gt.0) norm_ut = norm_ut + ghost_width_POST             ! everything's shifted by ghost_width_POST if not first process
            
            ! set plot variables
            grp_dim = [PB3D_plot%grid_eq%n(1),PB3D_plot%grid_eq%n(2),&
                &grid_X_trim%grp_n_r]
            tot_dim = [PB3D_plot%grid_eq%n(1),PB3D_plot%grid_eq%n(2),&
                &grid_X_trim%n(3)]
            grp_offset = [0,0,grid_X_trim%i_min-1]
            allocate(X_tot(grp_dim(1),grp_dim(2),grp_dim(3),6))
            allocate(Y_tot(grp_dim(1),grp_dim(2),grp_dim(3),6))
            allocate(Z_tot(grp_dim(1),grp_dim(2),grp_dim(3),6))
            do kd = 1,6
                X_tot(:,:,:,kd) = XYZ(:,:,norm_ut(1):norm_ut(2),1)
                Y_tot(:,:,:,kd) = XYZ(:,:,norm_ut(1):norm_ut(2),2)
                Z_tot(:,:,:,kd) = XYZ(:,:,norm_ut(1):norm_ut(2),3)
            end do
            allocate(var_names(2))
            
            ! calculate energy terms
            ierr = calc_E(PB3D_plot,.false.,X_id,E_pot,E_kin,E_pot_int,&
                &E_kin_int)
            CHCKERR('')
            
            ! point to the trimmed versions
            E_pot_trim => E_pot(:,:,norm_ut(1):norm_ut(2),:)
            E_kin_trim => E_kin(:,:,norm_ut(1):norm_ut(2),:)
            
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
                &grp_offset=[grp_offset,0],X=X_tot_trim,Y=Y_tot_trim,&
                &Z=Z_tot_trim,description=description)
            call plot_HDF5(var_names_kin,trim(file_name)//'_IM',&
                &imagpart(E_kin_trim),tot_dim=[tot_dim,2],&
                &grp_offset=[grp_offset,0],X=X_tot_trim,Y=Y_tot_trim,&
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
                &grp_offset=[grp_offset,0],X=X_tot_trim,Y=Y_tot_trim,&
                &Z=Z_tot_trim,description=description)
            call plot_HDF5(var_names_pot,trim(file_name)//'_IM',&
                &imagpart(E_pot_trim),tot_dim=[tot_dim,6],&
                &grp_offset=[grp_offset,0],X=X_tot_trim,Y=Y_tot_trim,&
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
                &[grp_dim,2])),tot_dim=[tot_dim,2],grp_offset=[grp_offset,0],&
                &X=X_tot_trim,Y=Y_tot_trim,Z=Z_tot_trim,description=description)
            call plot_HDF5(var_names,trim(file_name)//'_IM',imagpart(reshape(&
                &[sum(E_pot_trim(:,:,:,1:2),4),sum(E_pot_trim(:,:,:,3:6),4)],&
                &[grp_dim,2])),tot_dim=[tot_dim,2],grp_offset=[grp_offset,0],&
                &X=X_tot_trim,Y=Y_tot_trim,Z=Z_tot_trim,description=description)
            
            ! E
            var_names(1) = '1. potential energy'
            var_names(2) = '2. kinetic energy'
            file_name = trim(i2str(X_id))//'_E'
            description = 'total potential and kinetic energy'
            call plot_HDF5(var_names,trim(file_name)//'_RE',realpart(reshape(&
                &[sum(E_pot_trim,4),sum(E_kin_trim,4)],[grp_dim,2])),&
                &tot_dim=[tot_dim,2],grp_offset=[grp_offset,0],X=X_tot_trim,&
                &Y=Y_tot_trim,Z=Z_tot_trim,description=description)
            call plot_HDF5(var_names,trim(file_name)//'_IM',imagpart(reshape(&
                &[sum(E_pot_trim,4),sum(E_kin_trim,4)],[grp_dim,2])),&
                &tot_dim=[tot_dim,2],grp_offset=[grp_offset,0],X=X_tot_trim,&
                &Y=Y_tot_trim,Z=Z_tot_trim,description=description)
            
            call lvl_ud(-1)
            
            ! clean up
            nullify(E_kin_trim,E_pot_trim)
            nullify(X_tot_trim,Y_tot_trim,Z_tot_trim)
            call dealloc_grid(grid_X_trim)
        end if
    contains
        ! calculate the energy terms
        integer function calc_E(PB3D,B_aligned,X_id,E_pot,E_kin,E_pot_int,&
            &E_kin_int) result(ierr)
            use num_vars, only: use_pol_flux_F, grp_rank, glb_n_procs, &
                &ghost_width_POST
            use PB3D_vars, only: PB3D_type
            use eq_vars, only: vac_perm
            use utilities, only: c, con2dis
            use grid_ops, only: calc_int_vol, trim_grid
            use MPI_utilities, only: get_ser_var, wait_MPI, get_ghost_arr
            
            character(*), parameter :: rout_name = 'calc_E'
            
            ! input / output
            type(PB3D_type), intent(in) :: PB3D                                 ! PB3D variables
            logical, intent(in) :: B_aligned                                    ! whether grid is field-aligned
            integer, intent(in) :: X_id                                         ! nr. of Eigenvalue
            complex(dp), intent(inout), allocatable :: E_pot(:,:,:,:)           ! potential energy
            complex(dp), intent(inout), allocatable :: E_kin(:,:,:,:)           ! kinetic energy
            complex(dp), intent(inout), allocatable :: E_pot_int(:)             ! integrated potential energy
            complex(dp), intent(inout), allocatable :: E_kin_int(:)             ! integrated kinetic energy
            
            ! local variables
            integer :: kd                                                       ! counter
            integer :: i_lo, i_hi                                               ! upper and lower index for interpolation of eq grid to X grid
            integer :: grp_dim(3)                                               ! group dimension
            integer :: grp_dim_1(3)                                             ! group dimension with ghost region of width 1
            integer :: norm_ut(2)                                               ! untrimmed normal indices for trimmed grids
            type(grid_type) :: grid_X_trim                                      ! trimmed X grid
            real(dp) :: grp_r_eq                                                ! grp_r_F of X grid interpolated in eq grid
            real(dp), allocatable :: h22(:,:,:), g33(:,:,:), J(:,:,:)           ! interpolated h_FD(2,2), g_FD(3,3) and J_FD
            real(dp), allocatable :: kappa_n(:,:,:), kappa_g(:,:,:)             ! interpolated kappa_n and kappa_g
            real(dp), allocatable :: sigma(:,:,:)                               ! interpolated sigma
            real(dp), allocatable :: D2p(:,:,:), rho(:,:,:)                     ! interpolated Dpres_FD, rho
            real(dp), allocatable :: ang_1(:,:,:), ang_2(:,:,:), norm(:)        ! coordinates of grid in which to calculate total energy
            complex(dp), allocatable :: XUQ(:,:,:,:)                            ! X, U, Q_n and Q_g
            complex(dp), allocatable :: E_int_tot(:)                            ! integrated potential or kinetic energy of all processes
            
            ! set grp_dim
            grp_dim = [PB3D%grid_eq%n(1:2),PB3D%grid_X%grp_n_r]                 ! includes ghost regions of width ghost_width_POST
            
            ! allocate local variables
            allocate(h22(grp_dim(1),grp_dim(2),grp_dim(3)))
            allocate(g33(grp_dim(1),grp_dim(2),grp_dim(3)))
            allocate(J(grp_dim(1),grp_dim(2),grp_dim(3)))
            allocate(kappa_n(grp_dim(1),grp_dim(2),grp_dim(3)))
            allocate(kappa_g(grp_dim(1),grp_dim(2),grp_dim(3)))
            allocate(sigma(grp_dim(1),grp_dim(2),grp_dim(3)))
            allocate(D2p(grp_dim(1),grp_dim(2),grp_dim(3)))
            allocate(rho(grp_dim(1),grp_dim(2),grp_dim(3)))
            allocate(ang_1(grp_dim(1),grp_dim(2),grp_dim(3)))
            allocate(ang_2(grp_dim(1),grp_dim(2),grp_dim(3)))
            allocate(norm(grp_dim(3)))
            allocate(XUQ(grp_dim(1),grp_dim(2),grp_dim(3),4))
            allocate(E_pot(grp_dim(1),grp_dim(2),grp_dim(3),6))
            allocate(E_kin(grp_dim(1),grp_dim(2),grp_dim(3),2))
            allocate(E_pot_int(6))
            allocate(E_kin_int(2))
            
            ! iterate over all normal points in X grid and interpolate
            do kd = 1,grp_dim(3)
                ierr = con2dis(PB3D%grid_X%grp_r_F(kd),grp_r_eq,&
                    &PB3D%grid_eq%grp_r_F)
                CHCKERR('')
                i_lo = floor(grp_r_eq)
                i_hi = ceiling(grp_r_eq)
                
                h22(:,:,kd) = PB3D%met%h_FD(:,:,i_lo,c([2,2],.true.),0,0,0)+&
                    &(grp_r_eq-i_lo)*&
                    &(PB3D%met%h_FD(:,:,i_hi,c([2,2],.true.),0,0,0)&
                    &-PB3D%met%h_FD(:,:,i_lo,c([2,2],.true.),0,0,0))
                g33(:,:,kd) = PB3D%met%g_FD(:,:,i_lo,c([3,3],.true.),0,0,0)+&
                    &(grp_r_eq-i_lo)*(&
                    &PB3D%met%g_FD(:,:,i_hi,c([3,3],.true.),0,0,0)&
                    &-PB3D%met%g_FD(:,:,i_lo,c([3,3],.true.),0,0,0))
                J(:,:,kd) = PB3D%met%jac_FD(:,:,i_lo,0,0,0)+(grp_r_eq-i_lo)*&
                    &(PB3D%met%jac_FD(:,:,i_hi,0,0,0)&
                    &-PB3D%met%jac_FD(:,:,i_lo,0,0,0))
                kappa_n(:,:,kd) = PB3D%eq%kappa_n(:,:,i_lo)+(grp_r_eq-i_lo)*&
                    &(PB3D%eq%kappa_n(:,:,i_hi)-PB3D%eq%kappa_n(:,:,i_lo))
                kappa_g(:,:,kd) = PB3D%eq%kappa_g(:,:,i_lo)+(grp_r_eq-i_lo)*&
                    &(PB3D%eq%kappa_g(:,:,i_hi)-PB3D%eq%kappa_g(:,:,i_lo))
                sigma(:,:,kd) = PB3D%eq%sigma(:,:,i_lo)+(grp_r_eq-i_lo)*&
                    &(PB3D%eq%sigma(:,:,i_hi)-PB3D%eq%sigma(:,:,i_lo))
                D2p(:,:,kd) = PB3D%eq%pres_FD(i_lo,1)+(grp_r_eq-i_lo)*&
                    &(PB3D%eq%pres_FD(i_hi,1)-PB3D%eq%pres_FD(i_lo,1))
                rho(:,:,kd) = PB3D%eq%rho(i_lo)+(grp_r_eq-i_lo)*&
                    &(PB3D%eq%rho(i_hi)-PB3D%eq%rho(i_lo))
                if (B_aligned) then
                    if (use_pol_flux_F) then
                        ang_1(:,:,kd) = PB3D%grid_eq%theta_F(:,:,i_lo)+&
                            &(grp_r_eq-i_lo)*(PB3D%grid_eq%theta_F(:,:,i_hi)-&
                            &PB3D%grid_eq%theta_F(:,:,i_lo))                    ! theta
                    else
                        ang_1(:,:,kd) = PB3D%grid_eq%zeta_F(:,:,i_lo)+&
                            &(grp_r_eq-i_lo)*(PB3D%grid_eq%zeta_F(:,:,i_hi)-&
                            &PB3D%grid_eq%zeta_F(:,:,i_lo))                     ! zeta
                    end if
                    ang_2(:,:,kd) = PB3D%alpha                                  ! alpha
                else
                    ang_1(:,:,kd) = PB3D%grid_eq%theta_F(:,:,i_lo)+&
                        &(grp_r_eq-i_lo)*(PB3D%grid_eq%theta_F(:,:,i_hi)-&
                        &PB3D%grid_eq%theta_F(:,:,i_lo))                        ! theta
                    ang_2(:,:,kd) = PB3D%grid_eq%zeta_F(:,:,i_lo)+&
                        &(grp_r_eq-i_lo)*(PB3D%grid_eq%zeta_F(:,:,i_hi)-&
                        &PB3D%grid_eq%zeta_F(:,:,i_lo))                         ! zeta
                end if
                norm(kd) = PB3D%grid_eq%grp_r_F(i_lo)+(grp_r_eq-i_lo)*&
                    &(PB3D%grid_eq%grp_r_F(i_hi)-PB3D%grid_eq%grp_r_F(i_lo))
            end do
            
            ! calculate X, U, Q_n and Q_g
            do kd = 1,4
                ierr = calc_XUQ_ind(PB3D%grid_eq,PB3D%eq,PB3D%grid_X,PB3D%X,&
                    &X_id,kd,0._dp,XUQ(:,:,:,kd),met=PB3D%met)
                CHCKERR('')
            end do
            
            ! calc kinetic energy
            E_kin(:,:,:,1) = rho/h22*XUQ(:,:,:,1)*conjg(XUQ(:,:,:,1))
            E_kin(:,:,:,2) = rho*h22*J**2/g33*XUQ(:,:,:,2)*conjg(XUQ(:,:,:,2))
            
            ! calc potential energy
            E_pot(:,:,:,1) = 1._dp/vac_perm*1._dp/h22*XUQ(:,:,:,3)*&
                &conjg(XUQ(:,:,:,3))
            E_pot(:,:,:,2) = 1._dp/vac_perm*h22*J**2/g33*XUQ(:,:,:,4)*&
                &conjg(XUQ(:,:,:,4))
            E_pot(:,:,:,3) = -2*D2p*kappa_n*XUQ(:,:,:,1)*conjg(XUQ(:,:,:,1))
            E_pot(:,:,:,4) = -2*D2p*kappa_g*XUQ(:,:,:,1)*conjg(XUQ(:,:,:,2))
            E_pot(:,:,:,5) = -sigma*XUQ(:,:,:,4)*conjg(XUQ(:,:,:,1))
            E_pot(:,:,:,6) = sigma*XUQ(:,:,:,3)*conjg(XUQ(:,:,:,2))
            
            ! trim X grid
            ierr = trim_grid(PB3D%grid_X,grid_X_trim,&
                &shift_grid=ghost_width_POST)
            CHCKERR('')
            
            ! set up norm_ut including ghosted region of width 1
            norm_ut = [1,grid_X_trim%grp_n_r]
            if (grp_rank.gt.0) norm_ut = norm_ut + ghost_width_POST             ! everything's shifted by ghost_width_POST if not first process
            if (grp_rank.lt.glb_n_procs-1) norm_ut(2) = norm_ut(2)+1            ! ghost region of width 1
            
            ! set grp_dim ghosted with width 1
            grp_dim_1 = [PB3D%grid_eq%n(1:2),grid_X_trim%grp_n_r]               ! trimmed
            if (grp_rank.lt.glb_n_procs-1) grp_dim_1(3) = grp_dim_1(3)+1        ! ghost region of width 1 added (for integrals)
            
            ! integrate energy using ghosted variables
            ierr = calc_int_vol(ang_1(:,:,norm_ut(1):norm_ut(2)),&
                &ang_2(:,:,norm_ut(1):norm_ut(2)),norm(norm_ut(1):norm_ut(2)),&
                &J(:,:,norm_ut(1):norm_ut(2)),&
                &E_kin(:,:,norm_ut(1):norm_ut(2),:),E_kin_int)
            CHCKERR('')
            ierr = calc_int_vol(ang_1(:,:,norm_ut(1):norm_ut(2)),&
                &ang_2(:,:,norm_ut(1):norm_ut(2)),norm(norm_ut(1):norm_ut(2)),&
                &J(:,:,norm_ut(1):norm_ut(2)),&
                &E_pot(:,:,norm_ut(1):norm_ut(2),:),E_pot_int)
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
            E_kin = E_kin/sum(E_kin_int)
            E_pot = E_pot/sum(E_kin_int)
            E_kin_int = E_kin_int/sum(E_kin_int)
            E_pot_int = E_pot_int/sum(E_kin_int)
            
            ! deallocate variables
            call dealloc_grid(grid_X_trim)
        end function calc_E
    end function decompose_energy
end module sol_ops
