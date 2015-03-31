!------------------------------------------------------------------------------!
!   Operations on the equilibrium variables                                    !
!------------------------------------------------------------------------------!
module eq_ops
#include <PB3D_macros.h>
    use num_vars, only: pi, dp, max_str_ln
    use messages, only: print_ar_2, lvl_ud, writo
    use output_ops, only: print_GP_3D, draw_GP, &
        &print_GP_2D
    use str_ops, only: i2str, r2strt
    
    implicit none
    private
    public read_eq, calc_flux_q, prepare_RZL, calc_RZL, normalize_eq_vars, &
        &calc_norm_const, adapt_HEL_to_eq, calc_eq
    
    interface calc_RZL
        module procedure calc_RZL_ind, calc_RZL_arr
    end interface

contains
    ! reads the equilibrium input file
    integer function read_eq() result(ierr)
        use num_vars, only: eq_style, glb_rank, use_pol_flux_E
        use VMEC, only: read_VMEC
        use HELENA, only: read_HEL
        use grid_vars, only: n_r_eq
        
        character(*), parameter :: rout_name = 'read_eq'
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! only do this for the group master
        if (glb_rank.eq.0) then
            ! choose which equilibrium style is being used:
            !   1:  VMEC
            !   2:  HELENA
            select case (eq_style)
                case (1)                                                        ! VMEC
                    ierr = read_VMEC(n_r_eq,use_pol_flux_E)
                    CHCKERR('')
                case (2)                                                        ! HELENA
                    ierr = read_HEL(n_r_eq,use_pol_flux_E)
                    CHCKERR('')
                case default
                    err_msg = 'No equilibrium style associated with '//&
                        &trim(i2str(eq_style))
                    ierr = 1
                    CHCKERR(err_msg)
            end select
        end if
    end function read_eq
    
    ! calculate  the equilibrium  quantities on  a grid  determined by  straight
    ! field lines. This  grid has the dimensions  (n_par,grp_n_r_eq) where n_par
    ! is  the  number  of  points  taken along  the  magnetic  field  lines  and
    ! grp_n_r_eq .le.  n_r_eq is the  normal extent  in the equilibrium  grid of
    ! this rank. It is determined so  that the perturbation quantities that will
    ! be needed in this rank are fully covered, so no communication is necessary
    integer function calc_eq(grid_eq,eq,met,X,alpha) result(ierr)
        use eq_vars, only: dealloc_eq, eq_type
        use metric_ops, only: calc_g_C, calc_g_C, calc_T_VC, calc_g_V, &
            &calc_T_VF, calc_inv_met, calc_g_F, calc_jac_C, calc_jac_V, &
            &calc_jac_F, calc_f_deriv, calc_jac_H, calc_T_HF, calc_h_H
        use metric_vars, only: create_metric, dealloc_metric, metric_type
        use utilities, only: derivs
        use num_vars, only: max_deriv, ltest, plot_grid, eq_style
        use grid_vars, only: grid_type
        use grid_ops, only: calc_ang_grid, plot_grid_real
        use HELENA, only: dealloc_HEL
        use X_vars, only: X_type
        
        use utilities, only: calc_det, calc_inv, calc_mult, c
        !use metric_ops, only: plot_info
        !use utilities, only: plot_info_2 => plot_info
        
        character(*), parameter :: rout_name = 'calc_eq'
        
        ! input / output
        type(grid_type) :: grid_eq                                              ! equilibrium grid
        type(eq_type), intent(inout) :: eq                                      ! equilibrium variables
        type(metric_type), intent(inout) :: met                                 ! metric variables
        type(X_type), intent(inout) :: X                                        ! perturbation variables
        real(dp), intent(in) :: alpha                                           ! field line coordinate of current equilibrium
        
        ! local variables
        integer :: id
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer :: pmone                                                        ! plus or minus one
        
        !! TEMPORARILY, TO PLOT MORE INFO IN ROUTINES
        !plot_info = .true.
        !plot_info_2 = .true.
        
        !real(dp) :: A(1,1,1,16)
        !real(dp) :: B(1,1,1,10)
        !real(dp) :: det(1,1,1)
        !real(dp) :: invA(1,1,1,16)
        !real(dp) :: invB(1,1,1,10)
        !real(dp) :: AinvA(1,1,1,16)
        !real(dp) :: BinvB(1,1,1,10)
        !real(dp) :: AinvB(1,1,1,16)
        !real(dp) :: BinvA(1,1,1,16)
        
        !A = reshape([1,2,3,1,2,1,1,4,3,1,3,2,1,4,2,2],shape(A))
        !B = reshape([1,2,3,1,1,1,4,3,2,2],shape(B))
        !write(*,*) 'A = ', A
        !write(*,*) 'B = ', B
        !ierr = calc_det(det,A,4)
        !CHCKERR('')
        !write(*,*) 'detA = ', det
        !ierr = calc_det(det,B,4)
        !CHCKERR('')
        !write(*,*) 'detB = ', det
        !ierr = calc_inv(invA,A,4)
        !CHCKERR('')
        !write(*,*) 'invA = ', invA
        !ierr = calc_inv(invB,B,4)
        !CHCKERR('')
        !write(*,*) 'invB = ', invB
        !ierr = calc_mult(A,invA,AinvA,4)
        !CHCKERR('')
        !write(*,*) 'AinvA = ', AinvA
        !ierr = calc_mult(B,invB,BinvB,4)
        !CHCKERR('')
        !write(*,*) 'BinvB = ', BinvB
        !ierr = calc_mult(A,invB,AinvB,4)
        !CHCKERR('')
        !write(*,*) 'AinvB = ', AinvB
        !ierr = calc_mult(B,invA,BinvA,4)
        !CHCKERR('')
        !write(*,*) 'BinvA = ', BinvA
        
        ! initialize ierr
        ierr = 0
        
        call writo('Start setting up equilibrium quantities in '//&
            &trim(i2str(grid_eq%n(1)))//' discrete parallel points')
        call lvl_ud(1)
        
            call writo('Start initializing variables')
            call lvl_ud(1)
                
                ! initialize metric quantities
                call writo('Initialize metric quantities...')
                ierr = create_metric(grid_eq,met)
                CHCKERR('')
                
                ! calculate flux quantities and complete equilibrium grid
                call writo('Calculate flux quantities...')
                ierr = calc_flux_q(eq,grid_eq,X)
                CHCKERR('')
                
            call lvl_ud(-1)
            call writo('Variables initialized')
            
            call writo('Start determining the equilibrium grid')
            call lvl_ud(1)
            
                ! calculate  angular grid  points (theta,zeta)  that follow  the
                ! magnetic field  line determined  by alpha in  E(quilibrum) and
                ! F(lux) coords.
                ! Note: The  normal grid  is determined  by the  equilibrium, it
                ! can  either  use  the  toroidal  flux  or  the  poloidal  flux
                ! (use_pol_flux_E).  This  does  not   have  to  coincide  with
                ! use_pol_flux_F used by PB3D.
                ierr = calc_ang_grid(grid_eq,eq,alpha)
                CHCKERR('')
                
                ! adapt the HELENA variables to the equilibrium grid
                if (eq_style.eq.2) then
                    ierr = adapt_HEL_to_eq(grid_eq)
                    CHCKERR('')
                end if
                
                ! plot grid if requested
                if (plot_grid) then
                    ierr = plot_grid_real(grid_eq)
                    CHCKERR('')
                else
                    call writo('Grid plot not requested')
                end if
            
            call lvl_ud(-1)
            call writo('Equilibrium grid determined')
            
            call writo('Calculating equilibrium quantities on equilibrium grid')
            call lvl_ud(1)
            
                ! choose which equilibrium style is being used:
                !   1:  VMEC
                !   2:  HELENA
                select case (eq_style)
                    case (1)                                                    ! VMEC
                        ! calculate the  cylindrical variables  R, Z  and lambda
                        ! and derivatives
                        call writo('Calculate R,Z,L...')
                        ierr = prepare_RZL(grid_eq)
                        CHCKERR('')
                        do id = 0,max_deriv+1
                            ierr = calc_RZL(grid_eq,eq,derivs(id))
                            CHCKERR('')
                        end do
                        
                        ! calculate  the metrics  in the  cylindrical coordinate
                        ! system
                        call writo('Calculate g_C...')                          ! h_C is not necessary
                        do id = 0,max_deriv
                            ierr = calc_g_C(eq,met,derivs(id))
                            CHCKERR('')
                        end do
                        
                        ! calculate the  jacobian in the  cylindrical coordinate
                        ! system
                        call writo('Calculate jac_C...')
                        do id = 0,max_deriv
                            ierr = calc_jac_C(eq,met,derivs(id))
                            CHCKERR('')
                        end do
                        
                        ! calculate the  transformation matrix  C(ylindrical) ->
                        ! V(MEC)
                        call writo('Calculate T_VC...')
                        do id = 0,max_deriv
                            ierr = calc_T_VC(eq,met,derivs(id))
                            CHCKERR('')
                        end do
                        
                        ! calculate  the metric  factors in the  VMEC coordinate
                        ! system
                        call writo('Calculate g_V...')
                        do id = 0,max_deriv
                            ierr = calc_g_V(met,derivs(id))
                            CHCKERR('')
                        end do
                        
                        ! calculate the jacobian in the VMEC coordinate system
                        call writo('Calculate jac_V...')
                        do id = 0,max_deriv
                            ierr = calc_jac_V(met,derivs(id))
                            CHCKERR('')
                        end do
                        
                        ! calculate the transformation matrix V(MEC) -> F(lux)
                        call writo('Calculate T_VF...')
                        do id = 0,max_deriv
                            ierr = calc_T_VF(grid_eq,eq,met,derivs(id))
                            CHCKERR('')
                        end do
                        
                        ! set up plus minus one
                        pmone = -1                                              ! conversion VMEC LH -> RH coord. system
                    case (2)                                                    ! HELENA
                        ! calculate the jacobian in the HELENA coordinate system
                        call writo('Calculate jac_H...')
                        do id = 0,max_deriv
                            ierr = calc_jac_H(grid_eq,eq,met,derivs(id))
                            CHCKERR('')
                        end do
                        
                        ! calculate the metric factors  in the HELENA coordinate
                        ! system
                        call writo('Calculate h_H...')
                        do id = 0,max_deriv
                            ierr = calc_h_H(grid_eq,eq,met,derivs(id))
                            CHCKERR('')
                        end do
                        
                        ! calculate the inverse g_H of the metric factors h_H
                        call writo('Calculate g_H...')
                        do id = 0,max_deriv
                            ierr = calc_inv_met(met%g_E,met%h_E,derivs(id))
                            CHCKERR('')
                        end do
                        
                        ! calculate the transformation matrix H(ELENA) -> F(lux)
                        call writo('Calculate T_HF...')
                        do id = 0,max_deriv
                            ierr = calc_T_HF(grid_eq,eq,met,derivs(id))
                            CHCKERR('')
                        end do
                        
                        ! set up plus minus one
                        pmone = 1
                    case default
                        err_msg = 'No equilibrium style associated with '//&
                            &trim(i2str(eq_style))
                        ierr = 1
                        CHCKERR(err_msg)
                end select
                
                ! calculate  the  inverse of  the transformation  matrix T_EF
                call writo('Calculate T_FE...')
                do id = 0,max_deriv
                    ierr = calc_inv_met(met%T_FE,met%T_EF,derivs(id))
                    CHCKERR('')
                    ierr = calc_inv_met(met%det_T_FE,met%det_T_EF,derivs(id))
                    CHCKERR('')
                end do
                
                ! calculate the metric factors in the Flux coordinate system
                call writo('Calculate g_F...')
                do id = 0,max_deriv
                    ierr = calc_g_F(met,derivs(id))
                    CHCKERR('')
                end do
                
                ! calculate the inverse h_F of the metric factors g_F
                call writo('Calculate h_F...')
                do id = 0,max_deriv
                    ierr = calc_inv_met(met%h_F,met%g_F,derivs(id))
                    CHCKERR('')
                end do
                
                ! calculate the jacobian in the Flux coordinate system
                call writo('Calculate jac_F...')
                do id = 0,max_deriv
                    ierr = calc_jac_F(met,derivs(id))
                    CHCKERR('')
                end do
                
                ! calculate   the  derivatives  in  Flux  coordinates  from  the
                ! derivatives in VMEC coordinates
                call writo('Calculate derivatives in Flux coordinates...')
                do id = 0,max_deriv
                    ! g_FD
                    ierr = calc_f_deriv(met%g_F,met%T_FE,met%g_FD,max_deriv,&
                        &derivs(id))
                    CHCKERR('')
                    
                    ! h_FD
                    ierr = calc_f_deriv(met%h_F,met%T_FE,met%h_FD,max_deriv,&
                        &derivs(id))
                    CHCKERR('')
                    
                    ! jac_FD
                    ierr = calc_f_deriv(met%jac_F,met%T_FE,met%jac_FD,&
                        &max_deriv,derivs(id))
                    CHCKERR('')
                    
                    ! pres_FD
                    ierr = calc_f_deriv(eq%pres_E,&
                        &met%T_FE(1,1,:,c([2,1],.false.),:,0,0),&
                        &eq%pres_FD(:,id),max_deriv,id)
                    CHCKERR('')
                        
                    ! flux_p_FD
                    ierr = calc_f_deriv(eq%flux_p_E,&
                        &met%T_FE(1,1,:,c([2,1],.false.),:,0,0),&
                        &eq%flux_p_FD(:,id),max_deriv,id)
                    CHCKERR('')
                        
                    ! flux_t_FD
                    ierr = calc_f_deriv(eq%flux_t_E,&
                        &met%T_FE(1,1,:,c([2,1],.false.),:,0,0),&
                        &eq%flux_t_FD(:,id),max_deriv,id)
                    CHCKERR('')
                    eq%flux_t_FD(:,id) = pmone * eq%flux_t_FD(:,id)             ! multiply by plus minus one
                    
                    ! q_saf_FD
                    ierr = calc_f_deriv(eq%q_saf_E,&
                        &met%T_FE(1,1,:,c([2,1],.false.),:,0,0),&
                        &eq%q_saf_FD(:,id),max_deriv,id)
                    CHCKERR('')
                    eq%q_saf_FD(:,id) = pmone * eq%q_saf_FD(:,id)               ! multiply by plus minus one
                    
                    ! rot_t_FD
                    ierr = calc_f_deriv(eq%rot_t_E,&
                        &met%T_FE(1,1,:,c([2,1],.false.),:,0,0),&
                        &eq%rot_t_FD(:,id),max_deriv,id)
                    CHCKERR('')
                    eq%rot_t_FD(:,id) = pmone * eq%rot_t_FD(:,id)               ! multiply by plus minus one
                end do
                
                !write(*,*) 'det_T_FE:', 0,0,0
                !call print_HDF5('X','X',met%det_T_FE(:,:,:,0,0,0))
                !read(*,*)
                !do id = 1,9
                    !write(*,*) 'h_FD:', id
                    !call print_HDF5('X','X',met%h_FD(:,:,:,id,0,0,0))
                    !read(*,*)
                !end do
                
                ! deallocate unused equilibrium quantities
                if (.not.ltest) then
                    call writo('Deallocate unused equilibrium and metric &
                        &quantities...')
                    ierr = dealloc_metric(met)
                    CHCKERR('')
                    ! general equilibrium
                    ierr = dealloc_eq(eq)
                    CHCKERR('')
                    ! specific equilibrium
                    ! choose which equilibrium style is being used:
                    !   1:  VMEC
                    !   2:  HELENA
                    select case (eq_style)
                        case (1)                                                ! VMEC
                            ! nothing
                        case (2)                                                ! HELENA
                            call dealloc_HEL
                        case default
                            err_msg = 'No equilibrium style associated with '//&
                                &trim(i2str(eq_style))
                            ierr = 1
                            CHCKERR(err_msg)
                    end select
                end if
            
            call lvl_ud(-1)
            call writo('Equilibrium quantities calculated on equilibrium grid')
            
        call lvl_ud(-1)
        call writo('Done setting up equilibrium quantities')
    end function calc_eq
    
    ! prepare the cosine  and sine factors that are used  in the inverse Fourier
    ! transformation of R, Z and L and derivatives
    integer function prepare_RZL(grid) result(ierr)
        use fourier_ops, only: calc_trigon_factors
        use grid_vars, only: grid_type
        
        character(*), parameter :: rout_name = 'prepare_RZL'
        
        ! input / output
        type(grid_type), intent(inout) :: grid                                  ! grid for which to prepare the trigoniometric factors
        
        ! initialize ierr
        ierr = 0
        
        ! calculate trigoniometric factors using theta_E and zeta_E
        ierr = calc_trigon_factors(grid%theta_E,grid%zeta_E,grid%trigon_factors)
        CHCKERR('')
    end function prepare_RZL
    
    ! calculate  R, Z  and Lambda  and derivatives  in VMEC  coordinates at  the
    ! grid  points given  by  the  variables theta_E  and  zeta_E (contained  in
    ! trigon_factors) and at  every normal point. The  derivatives are indicated
    ! by the variable "deriv" which has 3 indices
    ! Note: There is no HELENA equivalent  because for HELENA simulations, R and
    ! Z are not necessary for calculation of the metric coefficients, and L does
    ! not exist.
    integer function calc_RZL_ind(grid,eq,deriv) result(ierr)
        use fourier_ops, only: fourier2real
        use VMEC, only: R_c, R_s, Z_c, Z_s, L_c, L_s
        use utilities, only: check_deriv
        use num_vars, only: max_deriv
        use eq_vars, only: eq_type
        use grid_vars, only: grid_type
        
        character(*), parameter :: rout_name = 'calc_RZL_ind'
        
        ! input / output
        type(grid_type), intent(in) :: grid                                     ! grid for which to prepare the trigoniometric factors
        type(eq_type), intent(inout) :: eq                                      ! equilibrium
        integer, intent(in) :: deriv(3)                                         ! derivatives
        
        ! initialize ierr
        ierr = 0
        
        ! check the derivatives requested
        ierr = check_deriv(deriv,max_deriv+1,'calc_RZL')
        CHCKERR('')
        
        ! calculate the variables R,Z and their angular derivative
        ierr = fourier2real(R_c(:,:,grid%i_min:grid%i_max,deriv(1)),&
            &R_s(:,:,grid%i_min:grid%i_max,deriv(1)),grid%trigon_factors,&
            &eq%R_E(:,:,:,deriv(1),deriv(2),deriv(3)),[deriv(2),deriv(3)])
        CHCKERR('')
        ierr = fourier2real(Z_c(:,:,grid%i_min:grid%i_max,deriv(1)),&
            &Z_s(:,:,grid%i_min:grid%i_max,deriv(1)),grid%trigon_factors,&
            &eq%Z_E(:,:,:,deriv(1),deriv(2),deriv(3)),[deriv(2),deriv(3)])
        CHCKERR('')
        ierr = fourier2real(L_c(:,:,grid%i_min:grid%i_max,deriv(1)),&
            &L_s(:,:,grid%i_min:grid%i_max,deriv(1)),grid%trigon_factors,&
            &eq%L_E(:,:,:,deriv(1),deriv(2),deriv(3)),[deriv(2),deriv(3)])
        CHCKERR('')
    end function calc_RZL_ind
    integer function calc_RZL_arr(grid,eq,deriv) result(ierr)
        use eq_vars, only: eq_type
        use grid_vars, only: grid_type
        
        character(*), parameter :: rout_name = 'calc_RZL_arr'
        
        ! input / output
        type(grid_type), intent(in) :: grid                                     ! grid for which to prepare the trigoniometric factors
        type(eq_type), intent(inout) :: eq                                      ! equilibrium
        integer, intent(in) :: deriv(:,:)
        
        ! local variables
        integer :: id
        
        do id = 1, size(deriv,2)
            ierr = calc_RZL_ind(grid,eq,deriv(:,id))
            CHCKERR('')
        end do
    end function calc_RZL_arr

    ! Adapt the HELENA  quantities for the equilibrium parallel  grid taking the
    ! correct poloidal HELENA range (and possibly multiples).
    ! Note: HELENA  equilibria are axisymmetric,  which eliminates the  need for
    ! zeta dependence. However, theta can depend  on angle_1 and angle_2 as well
    ! as on r (see discussion at the grid_type definition), the h_H_ij are 3D as
    ! well.
    integer function adapt_HEL_to_eq(grid) result(ierr)
        use num_vars, only: pi
        use HELENA, only: ias, chi_H, h_H_11_full, h_H_12_full, h_H_33_full, &
            &h_H_11, h_H_12, h_H_33
        use utilities, only: interp_fun
        use grid_vars, only: grid_type
        
        character(*), parameter :: rout_name = 'adapt_HEL_to_eq'
        
        ! input / output
        type(grid_type), intent(in) :: grid                                     ! grid to which to adapt the HELENA quantities
        
        ! local variables
        integer :: id, jd, kd                                                   ! counters
        real(dp) :: par_loc                                                     ! local parallel (= poloidal) point
        
        ! initialize ierr
        ierr = 0
        
        ! reallocate the arrays
        allocate(h_H_11(grid%n(1),grid%n(2),grid%grp_n_r))
        allocate(h_H_12(grid%n(1),grid%n(2),grid%grp_n_r))
        allocate(h_H_33(grid%n(1),grid%n(2),grid%grp_n_r))
        
        ! For every poloidal point, check  which half poloidal circle it belongs
        ! to.  If this  is a  bottom part  and HELENA  is symmetric  (ias =  0),
        ! the  quantities have  to  be taken  from  their symmetric  counterpart
        ! (2pi-theta) and  the metric factors  h_H_12 carry an  additional minus
        ! sign.
        ! loop over all angles of ang_2
        do jd = 1,grid%n(2)
            ! loop over all normal points of this rank
            do kd = 1,grid%grp_n_r
                ! loop over all angles of ang_1
                do id = 1,grid%n(1)
                    ! set the local poloidal point from theta
                    par_loc = grid%theta_E(id,jd,kd)
                    ! add or subtract  2pi to the parallel angle until  it is at
                    ! least 0 to get principal range 0..2pi
                    if (par_loc.lt.0._dp) then
                        do while (par_loc.lt.0._dp)
                            par_loc = par_loc + 2*pi
                        end do
                    else if (par_loc.gt.2*pi) then
                        do while (par_loc.gt.2._dp*pi)
                            par_loc = par_loc - 2*pi
                        end do
                    end if
                    ! Interpolate the  HELENA variables poloidally,  taking into
                    ! account the possible symmetry
                    if (ias.eq.0 .and. par_loc.gt.pi) then
                        ierr = interp_fun(h_H_11(id,jd,kd),&
                            &h_H_11_full(:,grid%i_min-1+kd),2*pi-par_loc,&
                            &x=chi_H)
                        CHCKERR('')
                        ierr = interp_fun(h_H_12(id,jd,kd),&
                            &-h_H_12_full(:,grid%i_min-1+kd),2*pi-par_loc,&
                            &x=chi_H)                                           ! change of sign
                        CHCKERR('')
                        ierr = interp_fun(h_H_33(id,jd,kd),&
                            &h_H_33_full(:,grid%i_min-1+kd),2*pi-par_loc,&
                            &x=chi_H)
                        CHCKERR('')
                    else
                        ierr = interp_fun(h_H_11(id,jd,kd),&
                            &h_H_11_full(:,grid%i_min-1+kd),par_loc,&
                            &x=chi_H)
                        CHCKERR('')
                        ierr = interp_fun(h_H_12(id,jd,kd),&
                            &h_H_12_full(:,grid%i_min-1+kd),par_loc,&
                            &x=chi_H)
                        CHCKERR('')
                        ierr = interp_fun(h_H_33(id,jd,kd),&
                            &h_H_33_full(:,grid%i_min-1+kd),par_loc,&
                            &x=chi_H)
                        CHCKERR('')
                    end if
                end do
            end do
        end do
    end function adapt_HEL_to_eq
    
    ! Calculates flux quantities  and normal derivatives in  the VMEC coordinate
    ! system. Also sets the normal coordinate in the equilibrium grid.
    integer function calc_flux_q(eq,grid_eq,X) result(ierr)
        use num_vars, only: eq_style, max_deriv, grp_nr, &
            &use_pol_flux_E, use_pol_flux_F, plot_flux_q
        use utilities, only: calc_deriv, calc_int
        use eq_vars, only: eq_type
        use grid_vars, only: grid_type
        use X_vars, only: X_type
        
        character(*), parameter :: rout_name = 'calc_flux_q'
        
        ! input / output
        type(eq_type), intent(inout) :: eq                                      ! equilibrium for this alpha
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid
        type(X_type), intent(inout) :: X                                        ! perturbation
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! choose which equilibrium style is being used:
        !   1:  VMEC
        !   2:  HELENA
        select case (eq_style)
            case (1)                                                            ! VMEC
                ierr = calc_flux_q_VMEC()
                CHCKERR('')
            case (2)                                                            ! HELENA
                ierr = calc_flux_q_HEL()
                CHCKERR('')
            case default
                err_msg = 'No equilibrium style associated with '//&
                    &trim(i2str(eq_style))
                ierr = 1
                CHCKERR(err_msg)
        end select
        
        ! plot flux quantities if requested
        call lvl_ud(1)
        if (plot_flux_q .and. grp_nr.eq.0) then                                 ! only first group because it is the same for all the groups
            ierr = flux_q_plot(eq,grid_eq)
            CHCKERR('')
        else
            call writo('Flux quantities plot not requested')
        end if
        call lvl_ud(-1)
    contains
        ! VMEC version
        integer function calc_flux_q_VMEC() result(ierr)
            use VMEC, only: iotaf, phi, phi_r, presf
            
            character(*), parameter :: rout_name = 'calc_flux_q_VMEC'
            
            ! local variables
            integer :: kd                                                       ! counter
            real(dp), allocatable :: Dflux_p_full(:)                            ! version of dflux_p/dp on full normal equilibrium grid (1..n(3))
            real(dp), allocatable :: flux_p_int_full(:)                         ! version of integrated flux_p on full normal grid (1..n(3))
            
            ! initialize ierr
            ierr = 0
            
            ! set up helper variables to calculate poloidal flux
            allocate(Dflux_p_full(grid_eq%n(3)),flux_p_int_full(grid_eq%n(3)))
            Dflux_p_full = iotaf*phi_r
            ierr = calc_int(Dflux_p_full,1.0_dp/(grid_eq%n(3)-1.0_dp),&
                &flux_p_int_full)
            CHCKERR('')
            
            ! poloidal flux: calculate using iotaf and phi, phi_r (easier to use
            ! full normal equilibrium grid flux_p because of the integral)
            eq%flux_p_E(:,1) = Dflux_p_full(grid_eq%i_min:grid_eq%i_max)
            eq%flux_p_E(:,0) = flux_p_int_full(grid_eq%i_min:grid_eq%i_max)
            do kd = 2,max_deriv+1
                ierr = calc_deriv(eq%flux_p_E(:,1),eq%flux_p_E(:,kd),&
                    &grid_eq%n(3)-1._dp,kd-1,1)
                CHCKERR('')
            end do
                
            ! toroidal flux: copy from VMEC and derive
            eq%flux_t_E(:,0) = phi(grid_eq%i_min:grid_eq%i_max)
            eq%flux_t_E(:,1) = phi_r(grid_eq%i_min:grid_eq%i_max)
            do kd = 2,max_deriv+1
                ierr = calc_deriv(eq%flux_t_E(:,1),eq%flux_t_E(:,kd),&
                    &grid_eq%n(3)-1._dp,kd-1,1)
                CHCKERR('')
            end do
            
            ! pressure: copy from VMEC and derive
            eq%pres_E(:,0) = presf(grid_eq%i_min:grid_eq%i_max)
            do kd = 1, max_deriv+1
                ierr = calc_deriv(eq%pres_E(:,0),eq%pres_E(:,kd),&
                    &grid_eq%n(3)-1._dp,kd,1)
                CHCKERR('')
            end do
            
            ! safety factor
            eq%q_saf_E(:,0) = 1.0_dp/iotaf(grid_eq%i_min:grid_eq%i_max)
            do kd = 1,max_deriv+1
                ierr = calc_deriv(eq%q_saf_E(:,0),eq%q_saf_E(:,kd),&
                    &grid_eq%n(3)-1._dp,kd,1)
                CHCKERR('')
            end do
            
            ! rot. transform
            eq%rot_t_E(:,0) = iotaf(grid_eq%i_min:grid_eq%i_max)
            do kd = 1,max_deriv+1
                ierr = calc_deriv(eq%rot_t_E(:,0),eq%rot_t_E(:,kd),&
                    &grid_eq%n(3)-1._dp,kd,1)
                CHCKERR('')
            end do
            
            ! max flux  of eq grid and of X grid and normal coord. of eq grid in
            ! Flux coordinates
            if (use_pol_flux_F) then
                X%max_flux_F = flux_p_int_full(grid_eq%n(3))
                eq%max_flux_F = X%max_flux_F
                grid_eq%r_F = flux_p_int_full/eq%max_flux_F
                grid_eq%grp_r_F = eq%flux_p_E(:,0)/eq%max_flux_F
            else
                X%max_flux_F = - phi(grid_eq%n(3))                              ! conversion VMEC LH -> RH coord. system
                eq%max_flux_F = X%max_flux_F
                grid_eq%r_F = - phi/eq%max_flux_F                               ! conversion VMEC LH -> RH coord. system
                grid_eq%grp_r_F = - eq%flux_t_E(:,0)/eq%max_flux_F              ! conversion VMEC LH -> RH coord. system
            end if
            
            ! max flux  of eq grid and of X grid and normal coord. of eq grid in
            ! Equilibrium coordinates
            if (use_pol_flux_E) then
                X%max_flux_E = flux_p_int_full(grid_eq%n(3))
                eq%max_flux_E = X%max_flux_E
                grid_eq%r_E = flux_p_int_full/eq%max_flux_E
                grid_eq%grp_r_E = eq%flux_p_E(:,0)/eq%max_flux_E
            else
                X%max_flux_E = phi(grid_eq%n(3))
                eq%max_flux_E = X%max_flux_E
                grid_eq%r_E = phi/eq%max_flux_E
                grid_eq%grp_r_E = eq%flux_t_E(:,0)/eq%max_flux_E
            end if
            
            ! deallocate helper variables
            deallocate(Dflux_p_full,flux_p_int_full)
        end function calc_flux_q_VMEC
        
        ! HELENA version
        integer function calc_flux_q_HEL() result(ierr)
            use HELENA, only: qs, flux_H, p0
            
            character(*), parameter :: rout_name = 'calc_flux_q_HEL'
            
            ! local variables
            integer :: kd                                                       ! counter
            real(dp), allocatable :: Dflux_t_full(:)                            ! version of dflux_t/dp on full normal equilibrium grid (1..n(3))
            real(dp), allocatable :: flux_t_int_full(:)                         ! version of integrated flux_t on full normal equilibrium grid (1..n(3))
            real(dp), allocatable :: flux_H_r(:)                                ! normal derivative of flux_H
            
            ! initialize ierr
            ierr = 0
            
            ! set up helper variables to calculate toroidal flux
            ! calculate normal derivative of flux_H
            allocate(flux_H_r(grid_eq%n(3)))
            ierr = calc_deriv(flux_H,flux_H_r,flux_H,1,1)
            CHCKERR('')
            allocate(Dflux_t_full(grid_eq%n(3)),flux_t_int_full(grid_eq%n(3)))
            Dflux_t_full = qs*flux_H_r
            ierr = calc_int(Dflux_t_full,flux_H,flux_t_int_full)
            CHCKERR('')
                
            ! poloidal flux: copy from HELENA and derive
            eq%flux_p_E(:,0) = flux_H(grid_eq%i_min:grid_eq%i_max)
            do kd = 1,max_deriv+1
                ierr = calc_deriv(eq%flux_p_E(:,0),eq%flux_p_E(:,kd),&
                    &eq%flux_p_E(:,0),kd,1)
                CHCKERR('')
            end do
            
            ! toroidal flux: calculate using  qs and flux_H, flux_H_r (easier to
            ! use full normal equilibrium grid flux_t because of the integral)
            eq%flux_t_E(:,1) = Dflux_t_full(grid_eq%i_min:grid_eq%i_max)
            eq%flux_t_E(:,0) = flux_t_int_full(grid_eq%i_min:grid_eq%i_max)
            do kd = 2,max_deriv+1
                ierr = calc_deriv(eq%flux_t_E(:,1),eq%flux_t_E(:,kd),&
                    &eq%flux_p_E(:,0),kd-1,1)
                CHCKERR('')
            end do
            
            ! pressure: copy from HELENA and derive
            eq%pres_E(:,0) = p0(grid_eq%i_min:grid_eq%i_max)
            do kd = 1, max_deriv+1
                ierr = calc_deriv(eq%pres_E(:,0),eq%pres_E(:,kd),&
                    &eq%flux_p_E(:,0),kd,1)
                CHCKERR('')
            end do
            
            ! safety factor
            eq%q_saf_E(:,0) = qs(grid_eq%i_min:grid_eq%i_max)
            do kd = 1,max_deriv+1
                ierr = calc_deriv(eq%q_saf_E(:,0),eq%q_saf_E(:,kd),&
                    &eq%flux_p_E(:,0),kd,1)
                CHCKERR('')
            end do
            
            ! rot. transform
            eq%rot_t_E(:,0) = 1.0_dp/qs(grid_eq%i_min:grid_eq%i_max)
            do kd = 1,max_deriv+1
                ierr = calc_deriv(eq%rot_t_E(:,0),eq%rot_t_E(:,kd),&
                    &eq%flux_p_E(:,0),kd,1)
                CHCKERR('')
            end do
            
            ! max flux  of eq grid and of X grid and normal coord. of eq grid in
            ! Flux coordinates
            if (use_pol_flux_F) then
                X%max_flux_F = flux_H(grid_eq%n(3))
                eq%max_flux_F = X%max_flux_F
                grid_eq%r_F = flux_H/eq%max_flux_F
                grid_eq%grp_r_F = eq%flux_p_E(:,0)/eq%max_flux_F
            else
                X%max_flux_F = flux_t_int_full(grid_eq%n(3))
                eq%max_flux_F = X%max_flux_F
                grid_eq%r_F = flux_t_int_full/eq%max_flux_F
                grid_eq%grp_r_F = eq%flux_t_E(:,0)/eq%max_flux_F
            end if
            
            ! max flux  of eq grid and of X grid and normal coord. of eq grid in 
            ! Equilibrium coordinates (uses poloidal flux by default)
            X%max_flux_E = flux_H(grid_eq%n(3))
            eq%max_flux_E = X%max_flux_E
            grid_eq%r_E = flux_H/eq%max_flux_E
            grid_eq%grp_r_E = eq%flux_p_E(:,0)/eq%max_flux_E
            
            ! deallocate helper variables
            deallocate(Dflux_t_full,flux_t_int_full)
        end function calc_flux_q_HEL
    end function calc_flux_q
    
    ! plots the flux quantities in the perturbation grid
    !   safety factor q_saf
    !   rotational transform rot
    !   pressure pres
    !   poloidal flux flux_p
    !   toroidal flux flux_t
    ! [MPI] Only first group
    integer function flux_q_plot(eq,grid_eq) result(ierr)
        use num_vars, only: eq_style, use_pol_flux_F, output_style, &
            &grp_rank, grp_n_procs, no_plots
        use eq_vars, only: eq_type
        use grid_vars, only: grid_type
        
        character(*), parameter :: rout_name = 'flux_q_plot'
        
        ! input / output
        type(eq_type), intent(in) :: eq                                         ! equilibrium for this alpha
        type(grid_type), intent(in) :: grid_eq                                  ! normal grid
        
        ! local variables
        integer :: id                                                           ! counter
        integer :: n_vars = 5                                                   ! nr. of variables to plot
        integer :: n_norm                                                       ! nr. of normal points in plot variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        real(dp), allocatable :: X_plot_2D(:,:)                                 ! x values of 2D plot
        real(dp), allocatable :: Y_plot_2D(:,:)                                 ! y values of 2D plot
        character(len=max_str_ln), allocatable :: plot_titles(:)                ! plot titles
        
        ! initialize ierr
        ierr = 0
        
        ! bypass plots if no_plots
        if (no_plots) return
        
        ! user output
        call writo('Plotting flux quantities')
        
        call lvl_ud(1)
        
        ! set up the number of normal grid points in plot variables
        n_norm = grid_eq%grp_n_r
        !if (output_style.eq.1) then
            !if (grp_rank.lt.grp_n_procs-1) n_norm = n_norm + 1                  ! extra ghost point
        !end if
        
        ! initialize X_plot_2D and Y_plot_2D
        allocate(X_plot_2D(n_norm,n_vars))
        allocate(Y_plot_2D(n_norm,n_vars))
        
        ! set up plot titles and file names
        allocate(plot_titles(n_vars))
        plot_titles(1) = 'safety factor []'
        plot_titles(2) = 'rotational transform []'
        plot_titles(3) = 'pressure [pa]'
        plot_titles(4) = 'poloidal flux [Tm^2]'
        plot_titles(5) = 'toroidal flux [Tm^2]'
        
        ! fill the 2D version of the plot
        ! choose which equilibrium style is being used:
        !   1:  VMEC
        !   2:  HELENA
        select case (eq_style)
            case (1)                                                            ! VMEC
                Y_plot_2D(:,1) = -eq%q_saf_E(:,0)                               ! conversion VMEC LH -> RH coord. system
                Y_plot_2D(:,2) = -eq%rot_t_E(:,0)                               ! conversion VMEC LH -> RH coord. system
                Y_plot_2D(:,3) = eq%pres_E(:,0)
                Y_plot_2D(:,4) = eq%flux_p_E(:,0)
                Y_plot_2D(:,5) = -eq%flux_t_E(:,0)                              ! conversion VMEC LH -> RH coord. system
            case (2)                                                            ! HELENA
                Y_plot_2D(:,1) = eq%q_saf_E(:,0)
                Y_plot_2D(:,2) = eq%rot_t_E(:,0)
                Y_plot_2D(:,3) = eq%pres_E(:,0)
                Y_plot_2D(:,4) = eq%flux_p_E(:,0)
                Y_plot_2D(:,5) = eq%flux_t_E(:,0)
            case default
                err_msg = 'No equilibrium style associated with '//&
                    &trim(i2str(eq_style))
                ierr = 1
                CHCKERR(err_msg)
        end select
        
        ! 2D normal variable (Y_plot_2D tabulated in eq. grid)
        if (use_pol_flux_F) then
            X_plot_2D(1:grid_eq%grp_n_r,1) = eq%flux_p_E(:,0)/eq%max_flux_E
        else
            X_plot_2D(1:grid_eq%grp_n_r,1) = eq%flux_t_E(:,0)/eq%max_flux_E
        end if
        do id = 2,n_vars
            X_plot_2D(:,id) = X_plot_2D(:,1)
        end do
        
        ! plot the output
        ! choose which output style is being used:
        !   1:  GNUPlot
        !   2:  HDF5
        select case (output_style)
            case (1)                                                            ! GNUPlot
                ierr = flux_q_plot_GP()
                CHCKERR('')
            case (2)                                                            ! HDF5
                ierr = flux_q_plot_HDF5()
                CHCKERR('')
            case default
                err_msg = 'No output style associated with '//&
                    &trim(i2str(output_style))
                ierr = 1
                CHCKERR(err_msg)
        end select
        
        ! deallocate
        deallocate(X_plot_2D,Y_plot_2D)
        
        call lvl_ud(-1)
    contains
        ! plots the pressure and fluxes in GNUplot
        integer function flux_q_plot_GP() result(ierr)
            use MPI_ops, only: wait_MPI
            use output_ops, only: merge_GP
            
            character(*), parameter :: rout_name = 'flux_q_plot_HDF5'
            
            ! local variables
            character(len=max_str_ln), allocatable :: file_name(:)              ! file_name
            character(len=max_str_ln), allocatable :: file_names(:,:)           ! file_names of all processes
            
            ! initialize ierr
            ierr = 0
            
            ! set up file names
            allocate(file_name(2))
            file_name(1) = 'pres'
            file_name(2) = 'flux'
            
            ! plot the individual 2D output of this process (except q_saf lotand
            ! prot_t, as they are already ted in plot_jq) ressure
            call writo('The safety factor and rotational transform are not &
                &plotted here. Instead, use the input variable "plot_jq".')
            call print_GP_2D(plot_titles(3),trim(file_name(1))//'_'//&
                &trim(i2str(grp_rank))//'.dat',Y_plot_2D(:,3),X_plot_2D(:,3),&
                &draw=.false.)
            ! fluxes
            call print_GP_2D(trim(plot_titles(4))//', '//trim(plot_titles(5)),&
                &trim(file_name(2))//'_'//trim(i2str(grp_rank))//'.dat',&
                &Y_plot_2D(:,4:5),X_plot_2D(:,4:5),draw=.false.)
            
            ! wait for all processes
            ierr = wait_MPI()
            CHCKERR('')
            
            ! draw plot together by group master
            if (grp_rank.eq.0) then
                ! set up file names in array
                allocate(file_names(grp_n_procs,2))
                do id = 1,grp_n_procs
                    file_names(id,1) = trim(file_name(1))//'_'//&
                        &trim(i2str(id-1))//'.dat'                              ! pressure
                    file_names(id,2) = trim(file_name(2))//'_'//&
                        &trim(i2str(id-1))//'.dat'                              ! fluxes
                end do
                ! merge files
                call merge_GP(file_names(:,1),trim(file_name(1))//'.dat',&
                    &delete=.true.)                                             ! pressure
                call merge_GP(file_names(:,2),trim(file_name(2))//'.dat',&
                    &delete=.true.)                                             ! fluxes
                ! draw
                call draw_GP(plot_titles(3),trim(file_name(1))//'.dat',1,.true.,&
                    &.false.)                                                   ! pressure
                call draw_GP(trim(plot_titles(4))//', '//trim(plot_titles(5)),&
                    &trim(file_name(2))//'.dat',2,.true.,.false.)               ! fluxes
            end if
        end function flux_q_plot_GP
        
        ! convert 2D plot to real plot in 3D and output in HDF5
        integer function flux_q_plot_HDF5() result(ierr)
            use output_ops, only: print_HDF5
            use num_vars, only: n_theta_plot, n_zeta_plot
            use grid_ops, only: calc_XYZ_grid
            use grid_vars, only: create_grid, destroy_grid
            use MPI_ops, only: get_ser_var
            
            character(*), parameter :: rout_name = 'flux_q_plot_HDF5'
            
            ! local variables
            integer :: kd                                                       ! counter
            real(dp), allocatable :: theta_plot(:,:,:), zeta_plot(:,:,:)        ! theta and zeta for 3D plot
            real(dp), allocatable :: X_plot_3D(:,:,:)                           ! x values of 3D plot
            real(dp), allocatable :: Y_plot_3D(:,:,:)                           ! y values of 3D plot
            real(dp), allocatable :: Z_plot_3D(:,:,:)                           ! z values of 3D plot
            real(dp), allocatable :: X_plot(:,:,:,:)                            ! x values of total plot
            real(dp), allocatable :: Y_plot(:,:,:,:)                            ! y values of total plot
            real(dp), allocatable :: Z_plot(:,:,:,:)                            ! z values of total plot
            real(dp), allocatable :: f_plot(:,:,:,:)                            ! values of variable of total plot
            integer :: plot_dim(4)                                              ! total plot dimensions
            integer :: plot_grp_dim(4)                                          ! group plot dimensions
            integer :: plot_offset(4)                                           ! plot offset
            type(grid_type) :: grid_plot                                        ! grid for plotting
            integer, allocatable :: tot_i_min(:)                                ! i_min of equilibrium grid of all processes
            character(len=max_str_ln) :: file_name                              ! file name
            
            ! initialize ierr
            ierr = 0
            
            ! set up file name
            file_name = 'flux_quantities'
            
            ! initialize theta_plot and zeta_plot
            allocate(theta_plot(n_theta_plot,n_zeta_plot,grid_eq%grp_n_r))
            if (n_theta_plot.eq.1) then
                theta_plot = 0.0_dp
            else
                do id = 1,n_theta_plot
                    theta_plot(id,:,:) = &
                        &pi+(id-1.0_dp)*2*pi/(n_theta_plot-1.0_dp)              ! starting from pi gives nicer plots
                end do
            end if
            ! zeta equidistant
            allocate(zeta_plot(n_theta_plot,n_zeta_plot,grid_eq%grp_n_r))
            if (n_zeta_plot.eq.1) then
                zeta_plot = 0.0_dp
            else
                do id = 1,n_zeta_plot
                    zeta_plot(:,id,:) = (id-1.0_dp)*2*pi/(n_zeta_plot-1.0_dp)
                end do
            end if
            
            ! create plot grid (only group quantities needed)
            ierr = create_grid(grid_plot,&
                &[n_theta_plot,n_zeta_plot,grid_eq%n(3)],&
                &[grid_eq%i_min,grid_eq%i_max])
            CHCKERR('')
            grid_plot%theta_E = theta_plot
            grid_plot%zeta_E = zeta_plot
            grid_plot%grp_r_E = grid_eq%grp_r_E
            
            ! get min_i of equilibrium grid
            ierr = get_ser_var([grid_eq%i_min],tot_i_min,scatter=.true.)
            CHCKERR('')
            
            ! set up plot_dim, plot_grp_dim and plot_offset
            plot_grp_dim = [n_theta_plot,n_zeta_plot,grid_plot%grp_n_r,n_vars]
            plot_dim = [n_theta_plot,n_zeta_plot,&
                &grid_plot%n(3)-tot_i_min(1)+1,n_vars]                          ! subtract i_min of first process
            plot_offset = [0,0,grid_plot%i_min-tot_i_min(1),n_vars]             ! subtract i_min of first process
            
            ! allocate 3D X, Y and Z
            allocate(X_plot_3D(plot_grp_dim(1),plot_grp_dim(2),plot_grp_dim(3)))
            allocate(Y_plot_3D(plot_grp_dim(1),plot_grp_dim(2),plot_grp_dim(3)))
            allocate(Z_plot_3D(plot_grp_dim(1),plot_grp_dim(2),plot_grp_dim(3)))
            
            ! calculate 3D X,Y and Z
            ierr = calc_XYZ_grid(grid_plot,X_plot_3D,Y_plot_3D,Z_plot_3D)
            CHCKERR('')
            
            ! set up total plot variables
            allocate(X_plot(plot_grp_dim(1),plot_grp_dim(2),plot_grp_dim(3),&
                &plot_grp_dim(4)))
            allocate(Y_plot(plot_grp_dim(1),plot_grp_dim(2),plot_grp_dim(3),&
                &plot_grp_dim(4)))
            allocate(Z_plot(plot_grp_dim(1),plot_grp_dim(2),plot_grp_dim(3),&
                &plot_grp_dim(4)))
            allocate(f_plot(plot_grp_dim(1),plot_grp_dim(2),plot_grp_dim(3),&
                &plot_grp_dim(4)))
            do id = 1,n_vars
                X_plot(:,:,:,id) = X_plot_3D
                Y_plot(:,:,:,id) = Y_plot_3D
                Z_plot(:,:,:,id) = Z_plot_3D
            end do
            do kd = 1,grid_plot%grp_n_r
                f_plot(:,:,kd,1) = Y_plot_2D(kd,1)                              ! safey factor
                f_plot(:,:,kd,2) = Y_plot_2D(kd,2)                              ! rotational transform
                f_plot(:,:,kd,3) = Y_plot_2D(kd,3)                              ! pressure
                f_plot(:,:,kd,4) = Y_plot_2D(kd,4)                              ! poloidal flux
                f_plot(:,:,kd,5) = Y_plot_2D(kd,5)                              ! toroidal flux
            end do
            
            ! print the output using HDF5
            call print_HDF5(plot_titles,file_name,f_plot,plot_dim,&
                &plot_grp_dim,plot_offset,X_plot,Y_plot,Z_plot,col=1,&
                &description='Flux quantities')
            
            ! deallocate and destroy grid
            deallocate(theta_plot,zeta_plot)
            deallocate(X_plot_3D,Y_plot_3D,Z_plot_3D)
            deallocate(X_plot,Y_plot,Z_plot,f_plot)
            call destroy_grid(grid_plot)
        end function flux_q_plot_HDF5
    end function flux_q_plot
    
    ! normalizes equilibrium quantities using the normalization constants
    subroutine normalize_eq_vars(eq)
        use eq_vars, only: eq_type, psi_0, pres_0
        
        ! local variables
        type(eq_type) :: eq                                                     ! equilibrium
        integer :: id                                                           ! counter
        
        ! scale the quantities
        eq%pres_FD = eq%pres_FD/pres_0
        eq%flux_p_FD = eq%flux_p_FD/psi_0
        eq%flux_t_FD = eq%flux_t_FD/psi_0
        eq%max_flux_E = eq%max_flux_E/psi_0
        eq%max_flux_F = eq%max_flux_F/psi_0
        
        ! scale  the  derivatives  by  psi_p_0
        do id = 1,size(eq%pres_FD,2)-1
            eq%pres_FD(:,id) = eq%pres_FD(:,id) * psi_0**(id)
            eq%flux_p_FD(:,id) = eq%flux_p_FD(:,id) * psi_0**(id)
            eq%flux_t_FD(:,id) = eq%flux_t_FD(:,id) * psi_0**(id)
            eq%q_saf_FD(:,id) = eq%q_saf_FD(:,id) * psi_0**(id)
            eq%rot_t_FD(:,id) = eq%rot_t_FD(:,id) * psi_0**(id)
        end do
    end subroutine normalize_eq_vars
    
    ! sets up normalization constants:
    !   R_0:    major radius (= average R on axis)
    !   rho_0:  mass density on axis (set up through input variable)
    !   pres_0: pressure on axis
    !   B_0:    average magnetic field (= sqrt(mu_0 pres_0))
    !   psi_0:  reference flux (= R_0^2 B_0)
    ! [MPI] only global master
    integer function calc_norm_const() result(ierr)
        use num_vars, only: glb_rank, eq_style, mu_0, use_normalization
        use eq_vars, only: T_0, B_0, pres_0, psi_0, R_0, rho_0
        
        character(*), parameter :: rout_name = 'calc_norm_const'
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        if (use_normalization) then
            ! user output
            call writo('Calculating the normalization constants')
            
            call lvl_ud(1)
            
            ! calculation
            if (glb_rank.eq.0) then
                ! choose which equilibrium style is being used:
                !   1:  VMEC
                !   2:  HELENA
                select case (eq_style)
                    case (1)                                                    ! VMEC
                        call calc_norm_const_VMEC
                    case (2)                                                    ! HELENA
                        call calc_norm_const_HEL
                    case default
                        err_msg = 'No equilibrium style associated with '//&
                            &trim(i2str(eq_style))
                        ierr = 1
                        CHCKERR(err_msg)
                end select
            end if
            
            ! Alfven velocity
            T_0 = sqrt(mu_0*rho_0)*R_0/B_0 
            
            call writo('Major radius    R_0 = '//trim(r2strt(R_0))//' m')
            call writo('Pressure        pres_0 = '//trim(r2strt(pres_0))//' Pa')
            call writo('Mass density    rho_0 = '//trim(r2strt(rho_0))&
                &//' kg/m^3')
            call writo('Magnetic field  B_0 = '//trim(r2strt(B_0))//' T')
            call writo('Magnetic flux   psi_0 = '//trim(r2strt(psi_0))//' Tm^2')
            call writo('Alfven time     T_0 = '//trim(r2strt(T_0))//' s')
            
            call lvl_ud(-1)
            
            ! user output
            call writo('Normalization constants calculated')
        else
            ! user output
            call writo('Normalization not used')
        end if
    contains 
        ! VMEC version
        subroutine calc_norm_const_VMEC
            use VMEC, only: R_c, presf
            
            ! set the major  radius as the average value of  R_E on the magnetic
            ! axis
            R_0 = R_c(0,0,1,0)
            
            ! rho_0 is set up through an input variable with the same name
            
            ! set pres_0 as pressure on axis
            pres_0 = presf(1)
            
            ! set the reference value for B_0 from B_0 = sqrt(mu_0 pres_0)
            B_0 = sqrt(pres_0 * mu_0)
            
            ! set reference flux
            psi_0 = R_0**2 * B_0
        end subroutine calc_norm_const_VMEC
        
        ! HELENA version
        subroutine calc_norm_const_HEL
            use HELENA, only: R_0_H, B_0_H
            
            ! set the major radius as the HELENA normalization parameter
            R_0 = R_0_H
            
            ! rho_0 is set up through an input variable with the same name
            
            !  set the  reference  value  for B_0  as  the HELENA  normalization
            ! parameter
            B_0 = B_0_H
            
            ! set pres_0 as B_0^2/mu_0
            pres_0 = B_0**2/mu_0
            
            ! set reference flux
            psi_0 = R_0**2 * B_0
        end subroutine calc_norm_const_HEL
    end function calc_norm_const
end module eq_ops
