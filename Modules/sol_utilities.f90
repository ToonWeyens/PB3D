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
    use grid_vars, only: grid_type
    use eq_vars, only: eq_1_type, eq_2_type
    use X_vars, only: X_1_type, modes_type
    use sol_vars, only: sol_type

    implicit none
    private
    public calc_XUQ, calc_tot_sol_vec, calc_loc_sol_vec, interp_V_spline
#if ldebug
    public debug_calc_XUQ_arr, debug_interp_V_spline
#endif
    
    ! global variables
#if ldebug
    logical :: debug_calc_XUQ_arr = .false.                                     !< plot debug information for calc_XUQ_arr \ldebug
    logical :: debug_interp_V_spline = .false.                                  !< debug information for interp_v_spline \ldebug
#endif
    
    ! interfaces
    
    !> \public Calculates the normal (\f$\cdot_n\f$) or geodesic (\f$\cdot_g\f$)
    !! components   of   the   plasma  perturbation   \f$\vec{\xi}\f$   or   the
    !! magnetic  perturbation \f$\vec{Q}  = \nabla  \times \left(\vec{x}  \times
    !! \vec{B}\right)\f$.
    !!
    !! Which variables are plot depends on \c XUQ_style \cite Weyens3D :
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
    !! The perturbation grid is assumed to  have the same angular coordinates as
    !! the equilibrium grid, and the normal coordinates correspond to either the
    !! equilibrium grid (\c X_grid_style 1),  the solution grid (\c X_grid_style
    !! 2) or the enriched equilibrium grid (\c X_grid_style 3).
    !!
    !! The output grid,  furthermore, has the angular part  corresponding to the
    !! equilibrium grid, and the normal part to the solution grid.
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
    integer function calc_XUQ_arr(mds_X,mds_sol,grid_eq,grid_X,grid_sol,eq_1,&
        &eq_2,X,sol,X_id,XUQ_style,time,XUQ,deriv) result(ierr)
        
        use num_vars, only: use_pol_flux_F, norm_disc_prec_sol, X_grid_style
        use num_utilities, only: con2dis, spline
        use X_utilities, only: trim_modes
#if ldebug
        use eq_vars, only: max_flux_F
        use num_utilities, only: calc_int
#endif
        
        character(*), parameter :: rout_name = 'calc_XUQ_arr'
        
        ! input / output
        type(modes_type), intent(in), target :: mds_X                           !< general modes variables in perturbation grid
        type(modes_type), intent(in), target :: mds_sol                         !< general modes variables in solution grid
        type(grid_type), intent(in) :: grid_eq                                  !< equilibrium grid
        type(grid_type), intent(in) :: grid_X                                   !< perturbation grid
        type(grid_type), intent(in) :: grid_sol                                 !< solution grid
        type(eq_1_type), intent(in) :: eq_1                                     !< flux equilibrium
        type(eq_2_type), intent(in) :: eq_2                                     !< metric equilibrium
        type(X_1_type), intent(in), target :: X                                 !< perturbation variables
        type(sol_type), intent(in) :: sol                                       !< solution variables
        integer, intent(in) :: X_id                                             !< nr. of Eigenvalue
        integer, intent(in) :: XUQ_style                                        !< whether to calculate \f$X\f$, \f$U\f$, \f$Q_n\f$ or \f$Q_g\f$
        real(dp), intent(in) :: time(:)                                         !< time range
        complex(dp), intent(inout) :: XUQ(:,:,:,:)                              !< \f$X\f$, \f$U\f$, \f$Q_n\f$ or \f$Q_g\f$
        logical, intent(in), optional :: deriv                                  !< return parallel derivative
        
        ! local variables
        integer :: m                                                            ! secondary mode number
        integer :: kd_tot                                                       ! kd in total grid
        integer :: ind                                                          ! index of mode in tables
        integer :: n_mod_tot                                                    ! total number of modes
        integer :: n_t                                                          ! number of time points
        integer :: n_loc, m_loc                                                 ! local mode numbers
        integer :: id, jd, kd                                                   ! counters
        integer :: XUQ_dims(4)                                                  ! dimensions of XUQ
        integer :: kdl_X(2), kdl_sol(2)                                         ! limits on normal index for a mode in X and sol grids
        integer :: id_lim_X(2), id_lim_sol(2)                                   ! limits on total modes
        integer, pointer :: sec_X_loc(:,:), sec_sol_loc(:,:)                    ! pointers to secondary mode variables
        real(dp), allocatable :: expon(:,:,:)                                   ! exponent in solution grid
        real(dp), allocatable :: par_fac(:)                                     ! multiplicative factor due to parallel derivative, without iu, in perturbation grid
        real(dp), allocatable :: jq(:)                                          ! iota or q, interpolated at perturbation grid
        real(dp), allocatable :: S(:,:,:), inv_J(:,:,:)                         ! S and 1/J, interpolated at perturbation grid
        real(dp), pointer :: r_X_loc(:), r_sol_loc(:)                           ! local r_X and r_sol for a mode
        real(dp), pointer :: theta_F(:,:,:), zeta_F(:,:,:)                      ! angles in solution normal grid
        complex(dp) :: sqrt_sol_val_norm                                        ! normalized sqrt(sol_val)
        complex(dp), allocatable :: fac_0(:,:,:), fac_1(:,:,:)                  ! factor to be multiplied with X and DX in perturbation grid
        complex(dp), allocatable :: XUQ_loc(:,:)                                ! complex XUQ without time at a normal point
        complex(dp), allocatable :: sol_vec(:,:)                                ! total solution vector in solution grid
        complex(dp), allocatable :: sol_vec_tot(:,:,:)                          ! total solution vector in solution grid and derivatives
        complex(dp), allocatable :: Dsol_vec(:,:)                               ! normal derivative of sol_vec in solution grid
        complex(dp), pointer :: U0(:,:,:), U1(:,:,:)                            ! U0 and U1 or parallel derivative
        logical :: extrap = .true.                                              ! whether extrapolation is used
        logical :: deriv_loc                                                    ! local copy of deriv
        character(len=max_str_ln) :: err_msg                                    ! error message
        
#if ldebug
        real(dp), allocatable :: sol_vec_ALT(:)                                 ! sol_vec calculated from Dsol_vec in solution grid
#endif
        
        ! initialize ierr
        ierr = 0
        
        ! set up dimensions
        n_t = size(time)
        XUQ_dims = shape(XUQ)
        
        ! set up local copy of deriv
        deriv_loc = .false.
        if (present(deriv)) deriv_loc = deriv
        
        ! tests
        if (grid_eq%n(1).ne.grid_X%n(1) .or. grid_eq%n(2).ne.grid_X%n(2)) then
            ierr = 1
            err_msg = 'Grids need to be compatible'
            CHCKERR(err_msg)
        end if
        if (grid_eq%n(1).ne.XUQ_dims(1) .or. grid_eq%n(2).ne.XUQ_dims(2) .or. &
            &grid_sol%loc_n_r.ne.XUQ_dims(3) .or. n_t.ne.XUQ_dims(4)) then
            ierr = 1
            err_msg = 'XUQ needs to have the correct dimensions'
            CHCKERR(err_msg)
        end if
        
        ! limit input modes to output modes
        ! (X grid should comprise sol grid)
        ierr = trim_modes(mds_X,mds_sol,id_lim_X,id_lim_sol)
        CHCKERR('')
        sec_X_loc => mds_X%sec(id_lim_X(1):id_lim_X(2),:)
        sec_sol_loc => mds_sol%sec(id_lim_sol(1):id_lim_sol(2),:)
        
        n_mod_tot = size(sec_X_loc,1)
        
        if (size(sec_X_loc,1).ne.size(sec_sol_loc,1)) then
            ierr = 1
            err_msg = 'modes X and sol have inconsistent sizes'
            CHCKERR(err_msg)
        end if
        
        ! set up normalized sqrt(sol_val)
        sqrt_sol_val_norm = sqrt(sol%val(X_id))
        if (ip(sqrt_sol_val_norm).gt.0) sqrt_sol_val_norm = - sqrt_sol_val_norm ! exploding solution, not the decaying one
        sqrt_sol_val_norm = sqrt_sol_val_norm / abs(sqrt_sol_val_norm)          ! normalize
        
        ! allocate helper variables
        allocate(jq(XUQ_dims(3)))
        allocate(S(XUQ_dims(1),XUQ_dims(2),XUQ_dims(3)))
        if ((XUQ_style.eq.3 .or. XUQ_style.eq.4)) &
            &allocate(inv_J(XUQ_dims(1),XUQ_dims(2),XUQ_dims(3)))
        allocate(XUQ_loc(XUQ_dims(1),XUQ_dims(2)))
        allocate(Dsol_vec(sol%n_mod,XUQ_dims(3)))
        
        ! interpolate eq->sol
        if (use_pol_flux_F) then
            ierr = spline(grid_eq%loc_r_F,eq_1%q_saf_FD(:,0),&
                &grid_sol%loc_r_F,jq,ord=norm_disc_prec_sol)
            CHCKERR('')
        else
            ierr = spline(grid_eq%loc_r_F,eq_1%rot_t_FD(:,0),&
                &grid_sol%loc_r_F,jq,ord=norm_disc_prec_sol)
            CHCKERR('')
        end if
        do jd = 1,XUQ_dims(2)
            do id = 1,XUQ_dims(1)
                ierr = spline(grid_eq%loc_r_F,eq_2%S(id,jd,:),&
                    &grid_sol%loc_r_F,S(id,jd,:),ord=norm_disc_prec_sol)
                CHCKERR('')
                if ((XUQ_style.eq.3 .or. XUQ_style.eq.4)) then
                    ierr = spline(grid_eq%loc_r_F,&
                        &1._dp/eq_2%jac_FD(id,jd,:,0,0,0),grid_sol%loc_r_F,&
                        &inv_J(id,jd,:),ord=norm_disc_prec_sol)
                    CHCKERR('')
                end if
            end do
        end do
        
        ! initialize XUQ
        XUQ = 0._dp
        
        ! set up normal derivative of X vec
        ! Note: In X_style 2 (fast), the mode  indices vary as a function as the
        ! normal coordinate.  All of  these contiguous  ranges are  tabulated in
        ! variable 'sec'.
        if (XUQ_style.eq.2 .or. XUQ_style.eq.4) then
            ! set up variables
            allocate(sol_vec_tot(n_mod_tot,grid_sol%loc_n_r,0:1))
            sol_vec_tot = 0._dp
            
            ! convert to total solution vector and derivate
            do id = 0,1
                ierr = calc_tot_sol_vec(mds_sol,grid_sol,sol%vec(:,:,X_id),&
                    &sol_vec,deriv=id)
                CHCKERR('')
                sol_vec_tot(:,:,id) = sol_vec
            end do
            
#if ldebug
            if (debug_calc_XUQ_arr) then
                do m = 1,n_mod_tot
                    call writo('For mode '//trim(i2str(m))//', testing &
                        &whether Dsol_vec is correct by comparing its integral &
                        &with original sol_vec')
                    kdl_sol = sec_sol_loc(m,2:3)
                    allocate(sol_vec_ALT(1:kdl_sol(2)-kdl_sol(1)+1))
                    sol_vec_ALT = 0._dp
                    ierr = calc_int(rp(sol_vec_tot(m,kdl_sol(1):kdl_sol(2),1)),&
                        &grid_sol%loc_r_F(kdl_sol(1):kdl_sol(2)),sol_vec_ALT)
                    CHCKERR('')
                    call print_ex_2D(['TEST_RE_Dsol_vec'],'',reshape(&
                        &[rp(sol_vec_tot(m,kdl_sol(1):kdl_sol(2),0)),&
                        &sol_vec_ALT],[kdl_sol(2)-kdl_sol(1)+1,2]),x=reshape(&
                        &grid_sol%loc_r_F(kdl_sol(1):kdl_sol(2)),&
                        &[kdl_sol(2)-kdl_sol(1)+1,1])/max_flux_F*2*pi)
                    ierr = calc_int(ip(sol_vec_tot(m,kdl_sol(1):kdl_sol(2),1)),&
                        &grid_sol%loc_r_F(kdl_sol(1):kdl_sol(2)),sol_vec_ALT)
                    CHCKERR('')
                    call print_ex_2D(['TEST_IM_Dsol_vec'],'',reshape(&
                        &[ip(sol_vec_tot(m,kdl_sol(1):kdl_sol(2),0)),&
                        &sol_vec_ALT],[kdl_sol(2)-kdl_sol(1)+1,2]),x=reshape(&
                        &grid_sol%loc_r_F(kdl_sol(1):kdl_sol(2)),&
                        &[kdl_sol(2)-kdl_sol(1)+1,1])/max_flux_F*2*pi)
                    deallocate(sol_vec_ALT)
                end do
            end if
#endif
            
            ! convert back to local solution vector
            ierr = calc_loc_sol_vec(mds_sol,grid_sol%i_min,sol_vec_tot(:,:,1),&
                &Dsol_vec)
            CHCKERR('')
            
            ! clean up
            deallocate(sol_vec,sol_vec_tot)
        else
            Dsol_vec = 0._dp
        end if
        
        ! interpolate or point angles
        select case (X_grid_style)
            case (1,3)                                                          ! equilibrium or enriched
                allocate(theta_F(XUQ_dims(1),XUQ_dims(2),XUQ_dims(3)))
                allocate(zeta_F(XUQ_dims(1),XUQ_dims(2),XUQ_dims(3)))
                do jd = 1,XUQ_dims(2)
                    do id = 1,XUQ_dims(1)
                        ierr = spline(grid_eq%loc_r_F,grid_eq%theta_F(id,jd,:),&
                            &grid_sol%loc_r_F,theta_F(id,jd,:))
                        CHCKERR('')
                        ierr = spline(grid_eq%loc_r_F,grid_eq%zeta_F(id,jd,:),&
                            &grid_sol%loc_r_F,zeta_F(id,jd,:))
                        CHCKERR('')
                    end do
                end do
            case (2)                                                            ! solution
                theta_F => grid_X%theta_F
                zeta_F => grid_X%zeta_F
        end select
        
        ! iterate over all secondary modes
        do m = 1,n_mod_tot
            ! test
            if (sec_X_loc(m,1).ne.sec_sol_loc(m,1)) then
                ierr = 1
                err_msg = 'For mode '//trim(i2str(m))//&
                    &', no consistency for modes in grid_X and grid_sol'
                CHCKERR(err_msg)
            end if
            
            ! set up and normal limits
            kdl_X = sec_X_loc(m,2:3)
            kdl_sol = sec_sol_loc(m,2:3)
            
            ! convert limits to local
            kdl_X = kdl_X - grid_X%i_min + 1
            kdl_sol = kdl_sol - grid_sol%i_min + 1
            
            ! get perturbation and solution normal ranges for this mode number
            r_X_loc => grid_X%loc_r_F(kdl_X(1):kdl_X(2))
            r_sol_loc => grid_sol%loc_r_F(kdl_sol(1):kdl_sol(2))
            
            ! set up table index
            ind = sec_sol_loc(m,4)
            
            ! Set up multiplicative factors fac_0 and fac_1 which are used in
            ! XUQ = fac_0*X + fac_1*DX with fac_0 and fac_1:
            !   1. fac_0 = 1 (i(nq-m) or i(n-iota m) for deriv)
            !      fac_1 = 0
            !   2. fac_0 = U_0 (DU_0 for deriv)
            !      fac_1 = U_1 (DU_1 for deriv)
            !   3. fac_0 = i(nq-m)/J or i(n-iota m)/J (no deriv)
            !      fac_1 = 0
            !   4. fac_0 = DU_0/J - S (no deriv)
            !      fac_1 = DU_1/J
            ! where par_fac = nq-m or n-iota m.
            ! These quantities are then multiplied by the exponent for each mode
            ! number.
            allocate(expon(XUQ_dims(1),XUQ_dims(2),kdl_sol(1):kdl_sol(2)))      ! normal subset of grid_sol%loc_r_F
            allocate(fac_0(XUQ_dims(1),XUQ_dims(2),kdl_sol(1):kdl_sol(2)))      ! normal subset of grid_sol%loc_r_F
            allocate(fac_1(XUQ_dims(1),XUQ_dims(2),kdl_sol(1):kdl_sol(2)))      ! normal subset of grid_sol%loc_r_F
            allocate(par_fac(kdl_sol(1):kdl_sol(2)))                            ! normal subset of grid_sol%loc_r_F
            
            ! calculate parallel multiplicative factor par_fac and exponent
            do kd = kdl_sol(1),kdl_sol(2)
                ! set up local mode numbers
                kd_tot = kd + grid_sol%i_min - 1
                n_loc = mds_sol%n(kd_tot,ind)
                m_loc = mds_sol%m(kd_tot,ind)
#if ldebug
                if (use_pol_flux_F .and. (m_loc.ne.sec_X_loc(m,1)) .or. &
                    &(.not.use_pol_flux_F .and. (n_loc.ne.sec_X_loc(m,1)))) then
                    ierr = 1
                    err_msg = 'mode numbers not consistent'
                    CHCKERR(err_msg)
                end if
#endif
                
                if (use_pol_flux_F) then
                    par_fac(kd) = n_loc*jq(kd) - m_loc
                else
                    par_fac(kd) = n_loc - m_loc*jq(kd)
                end if
                expon(:,:,kd) = n_loc*zeta_F(:,:,kd) - &
                    &m_loc*theta_F(:,:,kd)
            end do
            
            ! set up multiplicative factors in solution grid
            select case (XUQ_style)
                case (1)                                                        ! calculating X
                    if (deriv_loc) then                                         ! parallel derivative
                        do kd = kdl_sol(1),kdl_sol(2)
                            fac_0(:,:,kd) = iu*par_fac(kd)
                        end do
                    else
                        fac_0 = 1._dp
                    end if
                    fac_1 = 0._dp
                case (2)                                                        ! calculating U
                    ! interpolate X-> sol or point
                    select case (X_grid_style)
                        case (1,3)                                              ! equilibrium or enriched
                            allocate(U0(XUQ_dims(1),XUQ_dims(2),&
                                &kdl_sol(1):kdl_sol(2)))
                            allocate(U1(XUQ_dims(1),XUQ_dims(2),&
                                &kdl_sol(1):kdl_sol(2)))
                            if (deriv_loc) then                                 ! parallel derivative
                                ierr = interp_V_spline(&
                                    &X%DU_0(:,:,kdl_X(1):kdl_X(2),ind),U0,&
                                    &r_X_loc,r_sol_loc,extrap)
                                CHCKERR('')
                                ierr = interp_V_spline(&
                                    &X%DU_1(:,:,kdl_X(1):kdl_X(2),ind),U1,&
                                    &r_X_loc,r_sol_loc,extrap)
                                CHCKERR('')
                            else
                                ierr = interp_V_spline(&
                                    &X%U_0(:,:,kdl_X(1):kdl_X(2),ind),U0,&
                                    &r_X_loc,r_sol_loc,extrap)
                                CHCKERR('')
                                ierr = interp_V_spline(&
                                    &X%U_1(:,:,kdl_X(1):kdl_X(2),ind),U1,&
                                    &r_X_loc,r_sol_loc,extrap)
                                CHCKERR('')
                            end if
                        case (2)                                                ! solution
                            if (deriv_loc) then                                 ! parallel derivative
                                U0 => X%DU_0(:,:,kdl_X(1):kdl_X(2),ind)
                                U1 => X%DU_1(:,:,kdl_X(1):kdl_X(2),ind)
                            else
                                U0 => X%U_0(:,:,kdl_X(1):kdl_X(2),ind)
                                U1 => X%U_1(:,:,kdl_X(1):kdl_X(2),ind)
                            end if
                    end select
                    
                    ! assign
                    fac_0 = U0
                    fac_1 = U1
                    
                    ! clean up
                    select case (X_grid_style)
                        case (1,3)                                              ! equilibrium or enriched
                            deallocate(U0,U1)
                        case (2)                                                ! solution
                            ! do nothing
                    end select
                case (3)                                                        ! calculating Qn
                    if (deriv_loc) then                                         ! parallel derivative
                        ierr = 1
                        err_msg = 'For XUQ_style = '//trim(i2str(XUQ_style))//&
                            &' calculation of parallel derivative not supported'
                        CHCKERR('')
                    else
                        do kd = kdl_sol(1),kdl_sol(2)
                            fac_0(:,:,kd) = iu*par_fac(kd)*inv_J(:,:,kd)
                        end do
                    end if
                    fac_1 = 0._dp
                case (4)                                                        ! calculating Qg
                    ! interpolate X-> sol or point
                    err_msg = 'For XUQ_style = '//trim(i2str(XUQ_style))//&
                        &' calculation of parallel derivative not supported'
                    CHCKERR('')
                    select case (X_grid_style)
                        case (1,3)                                              ! equilibrium or enriched
                            allocate(U0(XUQ_dims(1),XUQ_dims(2),&
                                &kdl_sol(1):kdl_sol(2)))
                            allocate(U1(XUQ_dims(1),XUQ_dims(2),&
                                &kdl_sol(1):kdl_sol(2)))
                            if (deriv_loc) then                                 ! parallel derivative
                                ierr = 1
                                CHCKERR(err_msg)
                            else
                                ierr = interp_V_spline(&
                                    &X%DU_0(:,:,kdl_X(1):kdl_X(2),ind),U0,&
                                    &r_X_loc,r_sol_loc,extrap)
                                CHCKERR('')
                                ierr = interp_V_spline(&
                                    &X%DU_1(:,:,kdl_X(1):kdl_X(2),ind),U1,&
                                    &r_X_loc,r_sol_loc,extrap)
                                CHCKERR('')
                            end if
                        case (2)                                                ! solution
                            if (deriv_loc) then                                 ! parallel derivative
                                ierr = 1
                                CHCKERR(err_msg)
                            else
                                U0 => X%DU_0(:,:,kdl_X(1):kdl_X(2),ind)
                                U1 => X%DU_1(:,:,kdl_X(1):kdl_X(2),ind)
                            end if
                    end select
                    
                    ! assign
                    fac_0 = U0*inv_J(:,:,kdl_sol(1):kdl_sol(2)) - &
                        &S(:,:,kdl_sol(1):kdl_sol(2))
                    fac_1 = U1*inv_J(:,:,kdl_sol(1):kdl_sol(2))
                    
                    ! clean up
                    select case (X_grid_style)
                        case (1,3)                                              ! equilibrium or enriched
                            deallocate(U0,U1)
                        case (2)                                                ! solution
                            ! do nothing
                    end select
                case default
                    err_msg = 'No XUQ style associated with '//&
                        &trim(i2str(XUQ_style))
                    ierr = 1
                    CHCKERR(err_msg)
            end select
            
            ! iterate over all normal points in solution modes variables
            do kd = kdl_sol(1),kdl_sol(2)
                ! set up loc complex XUQ without time at this normal point
                XUQ_loc = exp(iu*expon(:,:,kd)) * &
                    &(sol%vec(ind,kd,X_id)*fac_0(:,:,kd) + &
                    &Dsol_vec(ind,kd)*fac_1(:,:,kd))
                
                ! iterate over time steps
                do id = 1,n_t
                    ! add current mode
                    XUQ(:,:,kd,id) = XUQ(:,:,kd,id) + XUQ_loc * &
                        &exp(iu*sqrt_sol_val_norm*time(id)*2*pi)
                end do
            end do
            
            ! clean up
            deallocate(fac_0,fac_1,expon,par_fac)
        end do
        
        ! clean up
        select case (X_grid_style)
            case (1,3)                                                          ! equilibrium or enriched
                deallocate(theta_F,zeta_F)
            case (2)                                                            ! solution
                ! do nothing
        end select
        nullify(r_X_loc,r_sol_loc)
        nullify(theta_F,zeta_F)
        nullify(U0,U1)
        nullify(sec_X_loc,sec_sol_loc)
    end function calc_XUQ_arr
    !> \private (time) individual version
    integer function calc_XUQ_ind(mds_X,mds_sol,grid_eq,grid_X,grid_sol,eq_1,&
        &eq_2,X,sol,X_id,XUQ_style,time,XUQ,deriv) result(ierr)
        
        character(*), parameter :: rout_name = 'calc_XUQ_ind'
        
        ! input / output
        type(modes_type), intent(in) :: mds_X                                   !< general modes variables in perturbation grid
        type(modes_type), intent(in) :: mds_sol                                 !< general modes variables in solution grid
        type(grid_type), intent(in) :: grid_eq                                  !< equilibrium grid
        type(grid_type), intent(in) :: grid_X                                   !< perturbation grid
        type(grid_type), intent(in) :: grid_sol                                 !< solution grid
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
        ierr = calc_XUQ_arr(mds_X,mds_sol,grid_eq,grid_X,grid_sol,eq_1,eq_2,X,&
            &sol,X_id,XUQ_style,[time],XUQ_arr,deriv=deriv)
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
    !! can be different from the local mode numbers for X style 2 (fast):
    !!  -  local: the  size of  the saved  perturbation and  solution variables,
    !!  prescribed  by the  user as  input:
    !!      - either directly through \c n_mod_X
    !!      -  or through  the limits  on the  modes using  \c min_sec_X  and \c
    !!      max_sec_X.
    !!  - total: all  the possible modes that can resonate  in the plasma, which
    !!  can be  different from the  local number for \c  X_style 2, as  the mode
    !!  numbers depend on the normal coordinate.
    !!
    !! If the output variable is not allocated, it is done here.
    !!
    !! Optionally  a  derivative  can  be requested,  depending  on  the  normal
    !! discretization order.
    !!
    !! \note Need to provide the grid in which mds is tabulated.
    !!
    !! \return ierr
    integer function calc_tot_sol_vec(mds,grid_sol,sol_vec_loc,sol_vec_tot,&
        &deriv) result(ierr)
        use num_vars, only: X_style, norm_disc_prec_sol
        use X_vars, only: n_mod_X
        use num_utilities, only: spline
        
        character(*), parameter :: rout_name = 'calc_tot_sol_vec'
        
        ! input / output
        type(modes_type), intent(in) :: mds                                     !< general modes variables
        type(grid_type), intent(in) :: grid_sol                                 !< solution grid
        complex(dp), intent(in) :: sol_vec_loc(:,:)                             !< solution vector for local resonating modes
        complex(dp), intent(inout), allocatable :: sol_vec_tot(:,:)             !< solution vector for all possible resonating modes
        integer, intent(in), optional :: deriv                                  !< order of derivative
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer :: n_mod_loc                                                    ! local number of modes that are used in the calculations
        integer :: n_mod_tot                                                    ! total number of modes
        integer :: ld                                                           ! counter
        integer :: kdl(2)                                                       ! kd limits
        integer :: ld_loc                                                       ! local ld
        integer :: n_r                                                          ! normal size of variables
        integer :: deriv_loc                                                    ! local derivative
        integer :: max_deriv                                                    ! maximum derivative
        integer :: min_range                                                    ! minimum range for derivative
        complex(dp), allocatable :: Dsol_vec_loc(:)                             ! local Dsol_vec_tot in solution grid
        
        ! initialize ierr
        ierr = 0
        
        ! set n_r
        n_r = size(sol_vec_loc,2)
        
        ! set local derivative
        deriv_loc = 0
        if (present(deriv)) deriv_loc = deriv
        
        ! test
        select case (norm_disc_prec_sol)
            case (1)                                                            ! linear
                max_deriv = 1
                min_range = 2
            case (2)                                                            ! Akita hermite
                max_deriv = 1
                min_range = 4
            case (3)                                                            ! spline
                max_deriv = 3
                min_range = 4
            case default
                ierr = 1
                err_msg = 'normal discretization precision for solution &
                    &cannot be higher than 3 or lower than 0'
                CHCKERR(err_msg)
        end select
        if (deriv_loc.lt.0 .or. deriv_loc.gt.max_deriv) then
            ierr = 1
            err_msg = 'deriv must be between 0 and '//trim(i2str(max_deriv))//&
                &' for normal discretization '//trim(i2str(norm_disc_prec_sol))
            CHCKERR(err_msg)
        endif
        
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
                
                ! copy or derivate
                if (deriv_loc.eq.0) then                                        ! copy
                    sol_vec_tot = sol_vec_loc
                else                                                            ! derivate
                    do ld = 1,size(sol_vec_tot,1)
                        ierr = spline(grid_sol%r_F,sol_vec_loc(ld,:),&
                            &grid_sol%r_F,sol_vec_tot(ld,:),&
                            &ord=norm_disc_prec_sol,deriv=deriv_loc)
                        CHCKERR('')
                    end do
                end if
            case (2)                                                            ! fast
                ! set local and total nr. of modes
                n_mod_loc = n_mod_X
                n_mod_tot = size(mds%sec,1)
                
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
                
                ! track all possible modes and copy or derivate
                sol_vec_tot = 0._dp
                do ld = 1,size(mds%sec,1)
                    ld_loc = mds%sec(ld,1)-minval(mds%sec(:,1),1)+1
                    kdl = [mds%sec(ld,2),mds%sec(ld,3)]-grid_sol%i_min+1        ! normal range in solution vector
                    kdl = max(1,min(kdl,n_r))                                   ! limit to sol_vec_loc range
                    if (deriv_loc.eq.0) then
                        sol_vec_tot(ld_loc,kdl(1):kdl(2)) = &
                            &sol_vec_loc(mds%sec(ld,4),kdl(1):kdl(2))
                    else
                        allocate(Dsol_vec_loc(kdl(2)-kdl(1)+1))
                        Dsol_vec_loc = 0._dp
                        
                        if (kdl(2)-kdl(1)+1.ge.min_range) then                  ! only calculate if enough points
                            ierr = spline(grid_sol%r_F(kdl(1):kdl(2)),&
                                &sol_vec_loc(mds%sec(ld,4),kdl(1):kdl(2)),&
                                &grid_sol%r_F(kdl(1):kdl(2)),Dsol_vec_loc,&
                                &ord=norm_disc_prec_sol,deriv=deriv_loc)
                            CHCKERR('')
                        end if
                        sol_vec_tot(ld_loc,kdl(1):kdl(2)) = Dsol_vec_loc
                        
                        deallocate(Dsol_vec_loc)
                    end if
                end do
        end select
    end function calc_tot_sol_vec
    
    !> Calculate  the  local version  of  the  solution  vector from  the  total
    !! version.
    !!
    !! \see See calc_tot_sol_vec().
    !!
    !! \return ierr
    integer function calc_loc_sol_vec(mds,i_min,sol_vec_tot,sol_vec_loc) &
        &result(ierr)
        use num_vars, only: X_style
        use X_vars, only: n_mod_X
        
        character(*), parameter :: rout_name = 'calc_loc_sol_vec'
        
        ! input / output
        type(modes_type), intent(in) :: mds                                     !< general modes variables
        integer, intent(in) :: i_min                                            !< \c i_min of grid in which variables are tabulated
        complex(dp), intent(in) :: sol_vec_tot(:,:)                             !< solution vector for all possible resonating modes
        complex(dp), intent(inout), allocatable :: sol_vec_loc(:,:)             !< solution vector for local resonating modes
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer :: n_mod_loc                                                    ! local number of modes that are used in the calculations
        integer :: n_mod_tot                                                    ! total number of modes
        integer :: ld                                                           ! counter
        integer :: kdl(2)                                                       ! kd limits
        integer :: ld_loc                                                       ! local ld
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
                n_mod_tot = size(mds%sec,1)
                
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
                do ld = 1,size(mds%sec,1)
                    ld_loc = mds%sec(ld,1)-minval(mds%sec(:,1),1)+1
                    kdl = [mds%sec(ld,2),mds%sec(ld,3)]-i_min+1                 ! range in sol_vec_loc
                    kdl = max(1,min(kdl,n_r))                                   ! limit to sol_vec_loc range
                    sol_vec_loc(mds%sec(ld,4),kdl(1):kdl(2)) = &
                        &sol_vec_tot(ld_loc,kdl(1):kdl(2))
                end do
        end select
    end function calc_loc_sol_vec
    
    !> Interpolation for a quantity V using splines.
    !!
    !! Optionally, a  variable 'was_copied'  is returned that  indicates whether
    !! the interpolation degraded to a plain copy.
    integer function interp_V_spline(V_i,V_o,r_i,r_o,extrap,was_copied) &
        &result(ierr)
        
        use num_utilities, only: spline
        
        character(*), parameter :: rout_name = 'interp_V_spline'
        
        ! input / output
        complex(dp), intent(in) :: V_i(:,:,:)                                   ! input variable
        complex(dp), intent(out) :: V_o(:,:,:)                                  ! output variable
        real(dp), intent(in) :: r_i(:)                                          ! input r
        real(dp), intent(in) :: r_o(:)                                          ! output r
        logical, intent(in) :: extrap                                           ! extrapolation possible
        logical, intent(out), optional :: was_copied                            ! inerpolation degraded to copy
        
        ! local variables
        integer :: id, jd, kd                                                   ! counters
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
#if ldebug
        ! test sizes
        if (size(r_i).ne.size(V_i,3)) then
            ierr = 1
            err_msg = 'r_i and V_i are not compatible'
            CHCKERR(err_msg)
        end if
        if (size(r_o).ne.size(V_o,3)) then
            ierr = 1
            err_msg = 'r_o and V_o are not compatible'
            CHCKERR(err_msg)
        end if
        if (size(V_i,1).ne.size(V_o,1) .or. size(V_i,2).ne.size(V_o,2)) then
            ierr = 1
            err_msg = 'V_i and V_o are not compatible'
            CHCKERR(err_msg)
        end if
#endif
        
        ! initialize was_copied
        if (present(was_copied)) was_copied = .false.
        
        ! interpolate depending on normal size
        select case (size(r_i))
            case (1)                                                            ! copy directly
                ! copy used
                if (present(was_copied)) was_copied = .true.
                do kd = 1,size(V_o,3)
                    V_o(:,:,kd) = V_i(:,:,1)
                end do
            case (2:3)                                                          ! linear
                do jd = 1,size(V_i,2)
                    do id = 1,size(V_i,1)
                        ierr = spline(r_i,V_i(id,jd,:),r_o,V_o(id,jd,:),&
                            &ord=1,extrap=extrap)
                        CHCKERR('')
                    end do
                end do
            case (4:)                                                           ! spline
                do jd = 1,size(V_i,2)
                    do id = 1,size(V_i,1)
                        ierr = spline(r_i,V_i(id,jd,:),r_o,V_o(id,jd,:),&
                            &ord=3,extrap=extrap)
                        CHCKERR('')
                    end do
                end do
#if ldebug
            case default
                ierr = 1
                CHCKERR('Need more than 0 points')
#endif
        end select
        
#if ldebug
        ! visualize
        do jd = 1,size(V_i,2)
            do id = 1,size(V_i,1)
                if (debug_interp_V_spline) then
                    write(*,*) '    id,jd', id, jd
                    call plot_comp('real part',r_i,rp(V_i(id,jd,:)),&
                        &r_o,rp(V_o(id,jd,:)))
                    call plot_comp('imaginary part',r_i,ip(V_i(id,jd,:)),&
                        &r_o,ip(V_o(id,jd,:)))
                end if
            end do
        end do
#endif
    contains
#if ldebug
        !> \private Plot to compare interpolation to original variables
        subroutine plot_comp(plot_name,x_i,y_i,x_o,y_o)
            ! input / output
            character(len=*), intent(in) :: plot_name                           ! name of plot
            real(dp), intent(in) :: x_i(:)                                      ! x input
            real(dp), intent(in) :: y_i(:)                                      ! y input
            real(dp), intent(in) :: x_o(:)                                      ! x output
            real(dp), intent(in) :: y_o(:)                                      ! y output
            
            ! local variables
            real(dp), allocatable :: x_plot(:,:)                                ! x for plot
            real(dp), allocatable :: y_plot(:,:)                                ! y for plot
            character(len=max_str_ln) :: var_names(2)                           ! variable names
            
            allocate(x_plot(max(size(x_i),size(x_o)),2))
            allocate(y_plot(max(size(x_i),size(x_o)),2))
            
            x_plot(1:size(x_i),1) = x_i
            x_plot(size(x_i):size(x_plot,1),1) = x_i(size(x_i))
            x_plot(1:size(x_o),2) = x_o
            x_plot(size(x_o):size(x_plot,1),2) = x_o(size(x_o))
            
            y_plot(1:size(x_i),1) = y_i(:)
            y_plot(size(x_i):size(y_plot,1),1) = y_i(size(x_i))
            y_plot(1:size(x_o),2) = y_o(:)
            y_plot(size(x_o):size(y_plot,1),2) = y_o(size(x_o))
            
            var_names(1) = 'input '//plot_name
            var_names(2) = 'output '//plot_name
            call print_ex_2D(var_names,'',y_plot,x=x_plot)
        end subroutine plot_comp
#endif
    end function interp_V_spline
end module sol_utilities
