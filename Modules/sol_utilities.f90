!------------------------------------------------------------------------------!
!> Numerical utilities related to the solution vectors.
!------------------------------------------------------------------------------!
module sol_utilities
#include <PB3D_macros.h>
#include <wrappers.h>
    use str_utilities
    use output_ops
    use messages
    use num_vars, only: dp, iu, max_str_ln, pi
    use grid_vars, only: grid_type, disc_type
    use eq_vars, only: eq_1_type, eq_2_type
    use X_vars, only: X_1_type
    use sol_vars, only: sol_type

    implicit none
    private
    public calc_XUQ, calc_tot_sol_vec, calc_loc_sol_vec
#if ldebug
    public debug_calc_XUQ_arr
#endif
    
    ! global variables
#if ldebug
    logical :: debug_calc_XUQ_arr = .false.                                     ! plot debug information for calc_XUQ_arr \ldebug
#endif
    
    ! interfaces
    
    !> \public Calculates the normal (\f$\cdot_n\f$) or geodesic (\f$\cdot_g\f$)
    !! components   of   the   plasma  perturbation   \f$\vec{\xi}\f$   or   the
    !! magnetic  perturbation \f$\vec{Q}  = \nabla  \times \left(\vec{x}  \times
    !! \vec{B}\right)\f$.
    !!
    !! Which variables are plot depends on \c XUQ_style \cite Weyens3D:
    !!  - \c XUQ_style = 1: \f$X\f$ (supports parallel derivative)
    !!  - \c XUQ_style = 2: \f$U\f$ (supports parallel derivative)
    !!  - \c XUQ_style = 3: \f$Q_n\f$
    !!  - \c XUQ_style = 4: \f$Q_g\f$
    !!
    !! \f$X\f$ and \f$U\f$ support the calculation  of the parallel derivative.
    !!
    !! For \f$Q_n\f$ and \f$Q_g\f$, the metric  variables have to be provided as
    !! well.
    !!
    !! As  this operation  requires  the equilibrium,  metric, perturbation  and
    !! solution  variables,  the  grid  in  which  this  happens  requires  some
    !! explanation:
    !!
    !!  - The equilibrium  and metric variables are tabulated  in an equilibrium
    !!  grid, which  needs to have the  same angular extent as  the perturbation
    !!  grid, in which the perturbation variables are tabulated.
    !!  - The normal extent of these two grids can be different, though, as well
    !!  as  the normal  extent  of  the solution  grid,  in  which the  solution
    !!  variables are tabulated.
    !!  - Currently,  the perturbation variables  are supposed to have  the same
    !!  angular extent  as the  equilibrium and metric  variables, and  the same
    !!  normal extent  as the  solution variables. The  output then  follows the
    !!  perturbation grid.
    !!
    !! \return ierr
    interface calc_XUQ
        !> \public
        module procedure calc_XUQ_ind
        !> \public
        module procedure calc_XUQ_arr
    end interface

contains
    !> \private (time) array version
    integer function calc_XUQ_arr(grid_eq,grid_X,eq_1,eq_2,X,sol,X_id,&
        &XUQ_style,time,XUQ,deriv) result(ierr)
        
        use num_vars, only: use_pol_flux_F, norm_disc_prec_sol
        use num_utilities, only: con2dis, spline3
        use grid_utilities, only: setup_interp_data, apply_disc
        use X_vars, only: sec_X_ind
#if ldebug
        use num_utilities, only: calc_int
#endif
        
        character(*), parameter :: rout_name = 'calc_XUQ_arr'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  !< equilibrium grid
        type(grid_type), intent(in) :: grid_X                                   !< perturbation grid
        type(eq_1_type), intent(in) :: eq_1                                     !< flux equilibrium
        type(eq_2_type), intent(in) :: eq_2                                     !< metric equilibrium
        type(X_1_type), intent(in) :: X                                         !< perturbation variables
        type(sol_type), intent(in) :: sol                                       !< solution variables
        integer, intent(in) :: X_id                                             !< nr. of Eigenvalue
        integer, intent(in) :: XUQ_style                                        !< whether to calculate \f$X\f$, \f$U\f$, \f$Q_n\f$ or \f$Q_g\f$
        real(dp), intent(in) :: time(:)                                         !< time range
        complex(dp), intent(inout) :: XUQ(:,:,:,:)                              !< \f$X\f$, \f$U\f$, \f$Q_n\f$ or \f$Q_g\f$
        logical, intent(in), optional :: deriv                                  !< return parallel derivative
        
        ! local variables
        integer :: n_t                                                          ! number of time points
        integer :: n_mod_loc                                                    ! local number of modes that are used in the calculations
        integer :: n_mod_tot                                                    ! total number of modes
        integer :: id, kd, ld                                                   ! counter
        character(len=max_str_ln) :: err_msg                                    ! error message
        complex(dp) :: sqrt_sol_val_norm                                        ! normalized sqrt(sol_val)
        logical :: deriv_loc                                                    ! local copy of deriv
        complex(dp), allocatable :: fac_0(:,:,:), fac_1(:,:,:)                  ! factor to be multiplied with X and DX
        complex(dp), allocatable :: XUQ_loc(:,:)                                ! complex XUQ without time at a normal point
        complex(dp), allocatable :: sol_vec_tot(:,:)                            ! total solution vector
        complex(dp), allocatable :: Dsol_vec(:,:)                               ! normal derivative of sol_vec
        complex(dp), allocatable :: Dsol_vec_tot(:,:)                           ! total normal derivative of sol_vec 
        complex(dp), allocatable :: Dsol_vec_loc(:)                             ! local Dsol_vec_tot
        real(dp), allocatable :: par_fac(:)                                     ! multiplicative factor due to parallel derivative, without iu
        real(dp), allocatable :: expon(:,:,:)                                   ! exponent
        real(dp), allocatable :: jq(:)                                          ! iota or q, interpolated at solution grid
        real(dp), allocatable :: S(:,:,:), inv_J(:,:,:)                         ! S and 1/J, interpolated at solution grid
        type(disc_type) :: norm_interp_data                                     ! data for normal interpolation
#if ldebug
        real(dp), allocatable :: sol_vec_ALT(:)                                 ! sol_vec calculated from Dsol_vec
#endif
        
        ! initialize ierr
        ierr = 0
        
        ! set up n_t
        n_t = size(time)
        
        ! set  up local and  total nr.  of modes, which  can be different  for X
        ! style 2 (see discussion in sol_utilities)
        n_mod_loc = sol%n_mod
        n_mod_tot = size(sec_X_ind,2)
        
        ! tests
        if (grid_eq%n(1).ne.grid_X%n(1) .or. grid_eq%n(2).ne.grid_X%n(2)) then
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
        
        ! set up local copy of deriv
        deriv_loc = .false.
        if (present(deriv)) deriv_loc = deriv
        
        ! set up normalized sqrt(sol_val)
        sqrt_sol_val_norm = sqrt(sol%val(X_id))
        if (ip(sqrt_sol_val_norm).gt.0) sqrt_sol_val_norm = - sqrt_sol_val_norm ! exploding solution, not the decaying one
        sqrt_sol_val_norm = sqrt_sol_val_norm / abs(sqrt_sol_val_norm)          ! normalize
        
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
        ! Note  that  if  the  resolution  for sol_vec  is  bad,  the  numerical
        ! derivative is  very inaccurate, therefore only  smooth (i.e. physical)
        ! solutions should be looked at.
        allocate(Dsol_vec(n_mod_loc,grid_X%loc_n_r))
        allocate(XUQ_loc(size(XUQ,1),size(XUQ,2)))
        allocate(fac_0(grid_X%n(1),grid_X%n(2),grid_X%loc_n_r))
        allocate(fac_1(grid_X%n(1),grid_X%n(2),grid_X%loc_n_r))
        allocate(par_fac(grid_X%loc_n_r))
        allocate(jq(grid_X%loc_n_r))
        allocate(S(grid_X%n(1),grid_X%n(2),grid_X%loc_n_r))
        if ((XUQ_style.eq.3 .or. XUQ_style.eq.4)) &
            &allocate(inv_J(grid_X%n(1),grid_X%n(2),grid_X%loc_n_r))
        allocate(expon(grid_X%n(1),grid_X%n(2),grid_X%loc_n_r))
        
        ! setup normal interpolation data for equilibrium grid
        ierr = setup_interp_data(grid_eq%loc_r_F,grid_X%loc_r_F,&
            &norm_interp_data,norm_disc_prec_sol)
        CHCKERR('')
        
        ! interpolate
        if (use_pol_flux_F) then
            ierr = apply_disc(eq_1%q_saf_FD(:,0),norm_interp_data,jq)
            CHCKERR('')
        else
            ierr = apply_disc(eq_1%rot_T_FD(:,0),norm_interp_data,jq)
            CHCKERR('')
        end if
        ierr = apply_disc(eq_2%S,norm_interp_data,S,3)
        CHCKERR('')
        if ((XUQ_style.eq.3 .or. XUQ_style.eq.4)) then
            ierr = apply_disc(1._dp/eq_2%jac_FD(:,:,:,0,0,0),norm_interp_data,&
                &inv_J,3)
            CHCKERR('')
        end if
        
        ! clean up normal interpolation data
        call norm_interp_data%dealloc()
        
        ! initialize XUQ
        XUQ = 0._dp
        
        ! set up normal derivative of X vec
        ! Note: The mode indices of the different modes can change as a function
        ! of the normal coordinate. Therefore, to correctly derive the same mode
        ! in  the normal  coordinate, this  mode has  to be  tracked along  this
        ! coordinate. 
        ! This is done by setting up  an total solution vector that contains all
        ! the possible  mode numbers n_mod_tot,  which is then  derived normally
        ! and finally the  corresponding local resonating modes  are copied back
        ! at each normal surface.
        if (XUQ_style.eq.2 .or. XUQ_style.eq.4) then
            ! set up variables
            allocate(sol_vec_tot(n_mod_tot,grid_X%loc_n_r))
            allocate(Dsol_vec_tot(n_mod_tot,grid_X%loc_n_r))
            allocate(Dsol_vec_loc(grid_X%loc_n_r))
            sol_vec_tot = 0._dp
            
            ! convert to total solution vector
            ierr = calc_tot_sol_vec(grid_X%i_min,sol%vec(:,:,X_id),sol_vec_tot)
            CHCKERR('')
            
            ! derive
            do ld = 1,n_mod_tot
                ierr = spline3(norm_disc_prec_sol,grid_X%loc_r_F,&
                    &sol_vec_tot(ld,:),grid_X%loc_r_F,dynew=Dsol_vec_loc)
                CHCKERR('')
                Dsol_vec_tot(ld,:) = Dsol_vec_loc
            end do
            
#if ldebug
            if (debug_calc_XUQ_arr) then
                do ld = 1,n_mod_tot
                    call writo('For mode '//trim(i2str(ld))//', testing &
                        &whether Dsol_vec is correct by comparing its integral &
                        &with original sol_vec')
                    allocate(sol_vec_ALT(1:grid_X%loc_n_r))
                    ierr = calc_int(rp(Dsol_vec_tot(ld,:)),&
                        &grid_X%loc_r_F,sol_vec_ALT)
                    CHCKERR('')
                    call print_ex_2D(['TEST_RE_Dsol_vec'],'',reshape(&
                        &[rp(sol_vec_tot(ld,:)),sol_vec_ALT],&
                        &[grid_X%loc_n_r,2]),x=reshape(&
                        &[grid_X%loc_r_F(1:grid_X%loc_n_r),&
                        &grid_X%loc_r_F(1:grid_X%loc_n_r)],&
                        &[grid_X%loc_n_r,2]))
                    ierr = calc_int(ip(Dsol_vec_tot(ld,:)),&
                        &grid_X%loc_r_F,sol_vec_ALT)
                    CHCKERR('')
                    call print_ex_2D(['TEST_IM_Dsol_vec'],'',reshape(&
                        &[ip(sol_vec_tot(ld,:)),sol_vec_ALT],&
                        &[grid_X%loc_n_r,2]),x=reshape(&
                        &[grid_X%loc_r_F(1:grid_X%loc_n_r),&
                        &grid_X%loc_r_F(1:grid_X%loc_n_r)],&
                        &[grid_X%loc_n_r,2]))
                    deallocate(sol_vec_ALT)
                end do
            end if
#endif
            
            ! convert back to local solution vector
            ierr = calc_loc_sol_vec(grid_X%i_min,Dsol_vec_tot,Dsol_vec)
            CHCKERR('')
            
            ! clean up
            deallocate(sol_vec_tot,Dsol_vec_tot)
        else
            Dsol_vec = 0._dp
        end if
        
        ! iterate over all Fourier modes
        Fourier: do ld = 1,n_mod_loc
            ! set up parallel multiplicative factor par_fac in normal eq grid
            if (use_pol_flux_F) then
                par_fac = sol%n(:,ld)*jq-sol%m(:,ld)
            else
                par_fac = sol%n(:,ld)-sol%m(:,ld)*jq
            end if
            
            ! set up exponent
            do kd = 1,grid_X%loc_n_r
                expon(:,:,kd) = sol%n(kd,ld)*grid_X%zeta_F(:,:,kd) - &
                    &sol%m(kd,ld)*grid_X%theta_F(:,:,kd)
            end do
            
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
                        fac_0 = X%DU_0(:,:,:,ld)
                        fac_1 = X%DU_1(:,:,:,ld)
                    else
                        fac_0 = X%U_0(:,:,:,ld)
                        fac_1 = X%U_1(:,:,:,ld)
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
                            fac_0(:,:,kd) = X%DU_0(:,:,kd,ld)*inv_J(:,:,kd) &
                                &- S(:,:,kd)
                            fac_1(:,:,kd) = X%DU_1(:,:,kd,ld)*inv_J(:,:,kd)
                        end do
                    end if
                case default
                    err_msg = 'No XUQ style associated with '//&
                        &trim(i2str(XUQ_style))
                    ierr = 1
                    CHCKERR(err_msg)
            end select
            
            ! iterate over all normal points in sol grid
            normal: do kd = 1,grid_X%loc_n_r
                ! set up loc complex XUQ without time at this normal point
                XUQ_loc = exp(iu*expon(:,:,kd)) * &
                    &(sol%vec(ld,kd,X_id)*fac_0(:,:,kd) + &
                    &Dsol_vec(ld,kd)*fac_1(:,:,kd))
                
                ! iterate over time steps
                do id = 1,n_t
                    ! add current mode
                    XUQ(:,:,kd,id) = XUQ(:,:,kd,id) + XUQ_loc * &
                        &exp(iu*sqrt_sol_val_norm*time(id)*2*pi)
                end do
            end do normal
        end do Fourier
    end function calc_XUQ_arr
    !> \private (time) individual version
    integer function calc_XUQ_ind(grid_eq,grid_X,eq_1,eq_2,X,sol,X_id,&
        &XUQ_style,time,XUQ,deriv) result(ierr)
        
        character(*), parameter :: rout_name = 'calc_XUQ_ind'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  !< equilibrium grid
        type(grid_type), intent(in) :: grid_X                                   !< perturbation grid
        type(eq_1_type), intent(in) :: eq_1                                     !< flux equilibrium
        type(eq_2_type), intent(in) :: eq_2                                     !< metric equilibrium
        type(X_1_type), intent(in) :: X                                         !< perturbation variables
        type(sol_type), intent(in) :: sol                                       !< solution variables
        integer, intent(in) :: X_id                                             !< nr. of Eigenvalue
        integer, intent(in) :: XUQ_style                                        !< whether to calculate X, U, Qn or Qg
        real(dp), intent(in) :: time                                            !< time
        complex(dp), intent(inout) :: XUQ(:,:,:)                                !< normal component of perturbation
        logical, intent(in), optional :: deriv                                  !< return parallel derivative
        
        ! local variables
        complex(dp), allocatable :: XUQ_arr(:,:,:,:)                            ! dummy variable to use array version
        
        ! allocate array version of XUQ
        allocate(XUQ_arr(size(XUQ,1),size(XUQ,2),size(XUQ,3),1))
        
        ! call array version
        ierr = calc_XUQ_arr(grid_eq,grid_X,eq_1,eq_2,X,sol,X_id,XUQ_style,&
            &[time],XUQ_arr,deriv=deriv)
        CHCKERR('')
        
        ! copy array to individual XUQ
        XUQ = XUQ_arr(:,:,:,1)
        
        ! deallocate array version of XUQ
        deallocate(XUQ_arr)
    end function calc_XUQ_ind
    
    !> Calculate  the  total version  of  the  solution  vector from  the  local
    !! version.
    !!
    !! This is the  solution vector for all of the  possible mode numbers, which
    !! can be different from the local mode numbers for X style 2:
    !!  -  local: the  size of  the saved  perturbation and  solution variables,
    !!  prescribed  by the  user as  input:
    !!  - either directly through \c n_mod_X
    !!      -  or through  the limits  on the  modes using  \c min_sec_X  and \c
    !!      max_sec_X.
    !!  - total: all  the possible modes that can resonate  in the plasma, which
    !!  can be  different from the  local number for \c  X_style 2, as  the mode
    !!  numbers depend on the normal coordinate.
    !!
    !! If the output variable is not allocated, it is done here.
    !!
    !! \note If the lowest limits of the  grid is not 1 (e.g. <tt>grid_X%i_min =
    !! 1</tt> for first  process), the input variable \c i_min  should be set to
    !! set correctly. For a full grid, it should be set to 1.
    !!
    !! \return ierr
    integer function calc_tot_sol_vec(i_min,sol_vec_loc,sol_vec_tot) &
        &result(ierr)
        use num_vars, only: X_style
        use X_vars, only: n_mod_X, sec_X_ind
        
        character(*), parameter :: rout_name = 'calc_tot_sol_vec'
        
        ! input / output
        integer, intent(in) :: i_min                                            !< \c i_min of grid in which variables are tabulated
        complex(dp), intent(in) :: sol_vec_loc(:,:)                             !< solution vector for local resonating modes
        complex(dp), intent(inout), allocatable :: sol_vec_tot(:,:)             !< solution vector for all possible resonating modes
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer :: n_mod_loc                                                    ! local number of modes that are used in the calculations
        integer :: n_mod_tot                                                    ! total number of modes
        integer :: kd, ld                                                       ! counters
        integer :: kd_loc                                                       ! local kd
        integer :: n_r                                                          ! normal size of variables
        
        ! set n_r
        n_r = size(sol_vec_loc,2)
        
        ! initialize ierr
        ierr = 0
        
        ! operations depending on X style
        select case (X_style)
            case (1)                                                            ! prescribed
                ! test allocation
                if (allocated(sol_vec_tot)) then
                    if (size(sol_vec_tot,1).ne.size(sol_vec_loc,1) .or. &
                        &size(sol_vec_tot,2).ne.size(sol_vec_loc,2)) then
                        ierr = 1
                        err_msg = 'Local and total solution vectors need to &
                            &be of the same length'
                        CHCKERR(err_msg)
                    end if
                else
                    allocate(sol_vec_tot(size(sol_vec_loc,1),n_r))
                end if
                sol_vec_tot = sol_vec_loc
            case (2)                                                            ! fast
                ! set local and total nr. of modes
                n_mod_loc = n_mod_X
                n_mod_tot = size(sec_X_ind,2)
                
                ! test allocation
                if (allocated(sol_vec_tot)) then
                    if (size(sol_vec_tot,1).ne.n_mod_tot .or. &
                        &size(sol_vec_tot,2).ne.size(sol_vec_loc,2)) then
                        ierr = 1
                        err_msg = 'Total solution vectors needs to be of &
                            &correct length '//trim(i2str(n_mod_tot))
                        CHCKERR(err_msg)
                    end if
                else
                    allocate(sol_vec_tot(n_mod_tot,n_r))
                end if
                
                ! track all possible modes and copy
                sol_vec_tot = 0._dp
                do ld = 1,n_mod_tot
                    do kd = 1,n_r
                        kd_loc = kd + i_min - 1
                        if (sec_X_ind(kd_loc,ld).ne.0) sol_vec_tot(ld,kd) = &
                            &sol_vec_loc(sec_X_ind(kd_loc,ld),kd)
                    end do
                end do
        end select
    end function calc_tot_sol_vec
    
    !> Calculate  the  local version  of  the  solution  vector from  the  total
    !! version.
    !!
    !! \see See calc_tot_sol_vec().
    !!
    !! \return ierr
    integer function calc_loc_sol_vec(i_min,sol_vec_tot,sol_vec_loc) &
        &result(ierr)
        use num_vars, only: X_style
        use X_vars, only: n_mod_X, sec_X_ind
        
        character(*), parameter :: rout_name = 'calc_loc_sol_vec'
        
        ! input / output
        integer, intent(in) :: i_min                                            !< \c i_min of grid in which variables are tabulated
        complex(dp), intent(in) :: sol_vec_tot(:,:)                             !< solution vector for all possible resonating modes
        complex(dp), intent(inout), allocatable :: sol_vec_loc(:,:)             !< solution vector for local resonating modes
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer :: n_mod_loc                                                    ! local number of modes that are used in the calculations
        integer :: n_mod_tot                                                    ! total number of modes
        integer :: kd, ld                                                       ! counters
        integer :: kd_loc                                                       ! local kd
        integer :: n_r                                                          ! normal size of variables
        
        ! set n_r
        n_r = size(sol_vec_tot,2)
        
        ! initialize ierr
        ierr = 0
        
        ! operations depending on X style
        select case (X_style)
            case (1)                                                            ! prescribed
                ! test allocation
                if (allocated(sol_vec_loc)) then
                    if (size(sol_vec_loc,1).ne.size(sol_vec_tot,1) .or. &
                        &size(sol_vec_loc,2).ne.size(sol_vec_tot,2)) then
                        ierr = 1
                        err_msg = 'Local and total solution vectors need to &
                            &be of the same length'
                        CHCKERR(err_msg)
                    end if
                else
                    allocate(sol_vec_loc(size(sol_vec_tot,1),n_r))
                end if
                sol_vec_loc = sol_vec_tot
            case (2)                                                            ! fast
                ! set local and total nr. of modes
                n_mod_loc = n_mod_X
                n_mod_tot = size(sec_X_ind,2)
                
                ! test allocation
                if (allocated(sol_vec_loc)) then
                    if (size(sol_vec_loc,1).ne.n_mod_loc .or. &
                        &size(sol_vec_loc,2).ne.size(sol_vec_tot,2)) then
                        ierr = 1
                        err_msg = 'Total solution vectors needs to be of &
                            &correct length '//trim(i2str(n_mod_loc))
                        CHCKERR(err_msg)
                    end if
                else
                    allocate(sol_vec_loc(n_mod_loc,n_r))
                end if
                
                ! track all possible modes and copy
                sol_vec_loc = 0._dp
                do kd = 1,n_r
                    kd_loc = kd + i_min - 1
                    do ld = 1,n_mod_tot
                        if (sec_X_ind(kd_loc,ld).ne.0) &
                            &sol_vec_loc(sec_X_ind(kd_loc,ld),kd) = &
                            &sol_vec_tot(ld,kd)
                    end do
                end do
        end select
    end function calc_loc_sol_vec
end module sol_utilities
